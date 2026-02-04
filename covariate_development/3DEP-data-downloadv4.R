# ============================================================
# Get 3DEP-data
#
# Dave White - 2/4/26
# ============================================================
# PURPOSE
#   This script builds a DEM that is *perfectly aligned* to a Mukey raster (from gNATSGO via soilDB),
#   while allowing the user to choose a final DEM cell size (final_res_m).
#   Alignment means: identical CRS, origin (grid anchor), and AOI-constrained extent, so that cell
#   boundaries coincide and any Mukey-derived analyses can be combined cell-by-cell without resampling.
#
# HIGH-LEVEL WORKFLOW
#   1) Read AOI, buffer it to reduce edge artifacts.
#   2) Pull Mukey tiles via WCS in a small grid to limit request size; build a UNION template grid.
#   3) Create a *final* Mukey-based raster grid at the requested resolution, snapped to Mukey origin.
#   4) Tile the buffered AOI; fetch 3DEP elevation tiles (via elevatr/AWS) in Mukey CRS; normalize NoData;
#      snap each DEM tile to the *final* grid (bilinear).
#   5) Stack the aligned DEM tiles and aggregate with NA-safe mean (per cell).
#   6) Write the final LZW-compressed GeoTIFF and run diagnostics on tile outputs.
#
# KEY DESIGN CHOICES
#   - Mukey alignment: categorical Mukey tiles are aligned via nearest neighbor; continuous DEM via bilinear.
#   - Final grid resolution is coerced to a multiple or divisor of the Mukey resolution to preserve alignment.
#   - Tile overlap is added to avoid edge resampling artifacts when snapping to the final grid.
#   - Robust fetching with retries to handle transient AWS/HTTP failures.
#   - Normalization guards against extreme Float32 sentinel values sometimes observed in cloud-hosted rasters.
#
# INPUTS (EDIT THESE IN THE "Paths & parameters" SECTION)
#   - tmp_dir:      local scratch directory for intermediate tiles
#   - aoi_path:     AOI polygon (Shapefile/GeoPackage/etc.)
#   - out_final_dem: path for final DEM GeoTIFF
#   - buffer_m:     AOI buffer distance in meters
#   - dem_z:        elevatr zoom level (approx. 5–10 m native source for 3DEP/Mapzen-style tiles)
#   - tile_n:       (ncol, nrow) tiling of the AOI to partition requests (larger grid => smaller requests)
#   - final_res_m:  desired final output raster resolution in meters
#
# OUTPUTS
#   - dem_tile_final_XX.tif: per-tile DEMs snapped to the *final* grid in tmp_dir
#   - dem_mukey_snapped_final.tif: final aggregated DEM aligned to Mukey CRS, origin, and AOI extent
#
# ASSUMPTIONS / REQUIREMENTS
#   - AOI geometry is valid or can be repaired with st_make_valid.
#   - soilDB::mukey.wcs is reachable (gNATSGO WCS) and returns tiles at 30 m (categorical Mukey).
#   - elevatr can fetch 3DEP via AWS with the specified CRS; GDAL/curl env vars are set for robust I/O.
#   - CRS is projected (meter units) to ensure buffer/overlap distances make sense; if AOI is geographic,
#     the code transforms to the Mukey CRS before DEM requests.
#
# PERFORMANCE NOTES
#   - ncpu_fetch uses (cores - 2) to stay responsive; adjust upward if tile downloads are slow.
#   - Increase tile_n if you observe timeouts; reduce if overhead is too high.
#   - dem_z=13 is a good balance for ~5–10 m source tiles; raising z increases download size and potential
#     for timeouts; lowering z increases resampling distance to reach final_res_m.
#
# REPRODUCIBILITY
#   - The script fixes alignment by using Mukey’s origin and extent as the anchor; grid math is deterministic.
#   - All resampling methods are explicitly specified (near/bilinear) to avoid defaults changing across versions.
#   - GDAL options are set for compression, tiling, and BigTIFF fallback.
#
# TROUBLESHOOTING
#   - If elevatr requests fail, raise retry_tries or retry_backoff_s, or increase tile_n to reduce request size.
#   - If you see extreme min/max in diagnostics (±1e20), confirm normalize_nodata is applied as written.
#   - If compareGeom fails at the end, check CRS/res/extent messages printed earlier to identify drift.
#
suppressPackageStartupMessages({
  pkgs <- c("elevatr","terra","sf","soilDB","purrr")
  new  <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
  if (length(new)) install.packages(new)
  lapply(pkgs, require, character.only = TRUE)
})

# ---------------------------
# Paths & parameters (EDIT THESE)
# ---------------------------
# NOTE:
#   - tmp_dir should be a fast local disk (e.g., SSD) with enough space for temporary tiles.
#   - aoi_path can be any sf-readable vector dataset; CRS will be validated/transformed later.
#   - out_final_dem is the final aligned product; overwritten if it exists (see writeRaster wopt).
#   - buffer_m helps reduce resampling edge artifacts when clipping/fetching; 500 m is conservative.
#   - dem_z determines native tile resolution fetched by elevatr; z=13 ~ 5–10 m tiles in many areas.
#     https://github.com/tilezen/joerd/blob/master/docs/data-sources.md#what-is-the-ground-resolution
#     z=15 ~ 5m, z=14 ~ 10m, 
#   - tile_n controls fetch granularity; increase for smaller per-request footprint (helps avoid timeouts).
#   - final_res_m should ideally be a multiple or divisor of Mukey res (often 30 m) to preserve alignment.
tmp_dir         <- "D:/federal_lands/26_Spodic_Intensity/temp"
aoi_path        <- "D:/federal_lands/26_Spodic_Intensity/SpodicIntensityDSM_Boundary/SpodosolHuc12Opt2.shp"
out_final_dem   <- "D:/federal_lands/26_Spodic_Intensity/dem_mukey_snapped_final.tif"

buffer_m        <- 500       # AOI buffer (meters) — reduces edge effects and ensures coverage beyond AOI boundary
dem_z           <- 13        # elevatr zoom (~5–10 m source tiles) — higher z => finer native tiles, larger downloads
tile_n          <- c(3, 3)   # AOI tile grid (increase if timeouts) — partitions AOI into 9 tile requests
retry_tries     <- 3         # number of retry attempts per tile fetch
retry_backoff_s <- 4         # seconds to wait between retries; exponential/backoff could be added if needed
ncpu_fetch      <- max(1, parallel::detectCores() - 2)  # parallel workers used by elevatr

# >>> Desired final output resolution (meters) <<<
# IMPORTANT:
#   - final_res_m is coerced to a clean multiple/divisor of Mukey resolution if needed to maintain exact alignment.
#   - Example good choices: 10, 30, 60, 90 (commonly multiples/divisors of 30 m Mukey).
final_res_m     <- 30        # e.g., 10, 30, 60, 90 (prefer multiples/divisors of Mukey res)

dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------
# GDAL / curl environment
# ---------------------------
# These environment variables optimize remote reads via VSI (virtual file systems) and HTTP:
#   - GDAL_DISABLE_READDIR_ON_OPEN: avoids expensive directory listing for remote sources
#   - CPL_VSIL_CURL_ALLOWED_EXTENSIONS: restricts to GeoTIFF files over HTTP/S3 for performance
#   - CPL_VSIL_CURL_NON_PREFETCH: turn off prefetch to reduce unnecessary reads
#   - GDAL_HTTP_TIMEOUT: fail more predictably under slow connections
#   - AWS_NO_SIGN_REQUEST: enables anonymous access to public AWS 3DEP buckets
Sys.setenv(
  GDAL_DISABLE_READDIR_ON_OPEN = "TRUE",
  CPL_VSIL_CURL_ALLOWED_EXTENSIONS = ".tif,.tiff",
  CPL_VSIL_CURL_NON_PREFETCH = "YES",
  GDAL_HTTP_TIMEOUT = "600",
  AWS_NO_SIGN_REQUEST = "YES"
)

options(timeout = 600)  # general R download timeout safeguard

terraOptions(progress = 1)  # show progress bars for terra operations (use 0 to silence)

# ---------------------------
# Helpers
# ---------------------------

# Map extreme Float32 sentinels to NA, then clamp to plausible NM elevations.
# RATIONALE:
#   - Some cloud-hosted rasters contain extreme sentinel values (e.g., ~±1e29) in Float32 tiles.
#   - First line converts those extremes to NA; second line optionally clamps to plausible elevation bounds,
#     preventing skew in summary stats and downstream models.
normalize_nodata <- function(r) {
  r <- ifel(r < -1e20 | r >  1e20, NA, r)   # catches ~±1e29 values you saw
  r <- ifel(r < -500  | r >  6000, NA, r)   # optional plausibility clamp; adjust if your AOI warrants
  r
}

# Align categorical Mukey tiles (nearest neighbor)
# WHY "near"? Mukey is categorical; bilinear would create unphysical blended categories.
align_mukey_cat <- function(x, target) terra::project(x, target, method = "near")

# Align DEM tiles (bilinear)
# WHY "bilinear"? DEM is continuous; bilinear preserves smoothness when resampling to the final grid.
align_dem_bilinear <- function(x, target) terra::project(x, target, method = "bilinear")

# NA-safe mean reducer for stacking
# Used in terra::app to combine multiple aligned DEM tiles per cell.
mean_na <- function(x) mean(x, na.rm = TRUE)

# Build a final grid from a Mukey raster via factor-preserving coarsen/disaggregate.
# This guarantees perfect alignment to Mukey’s origin & extent.
# PROCESS:
#   - Read Mukey resolution/CRS/extent.
#   - Adjust desired_res_m to nearest multiple/divisor if needed (exact grid alignment).
#   - Compute ncols/nrows so AOI extent is fully covered (ceil ensures xmax/ymax coverage).
#   - Create an empty raster with NA values; later DEM tiles are snapped onto this grid.
make_final_grid_from_mukey <- function(mukey_union, desired_res_m) {
  mu_res <- terra::res(mukey_union)[1]
  mu_crs <- terra::crs(mukey_union)
  mu_ext <- terra::ext(mukey_union)
  eps    <- 1e-6
  
  # Snap desired_res_m to a clean multiple/divisor of Mukey res
  ratio <- desired_res_m / mu_res
  if (abs(ratio - round(ratio)) > eps && abs(1/ratio - round(1/ratio)) > eps) {
    # pick nearest multiple or divisor to preserve alignment
    m_up <- round(ratio); m_down <- round(1/ratio)
    if (abs(ratio - m_up) <= abs(1/ratio - m_down)) desired_res_m <- m_up * mu_res
    else desired_res_m <- mu_res / m_down
    message(sprintf("Adjusted final resolution to %.3f m for exact Mukey alignment.", desired_res_m))
  }
  
  # Compute ncols/nrows by rounding UP so xmax/ymax are fully covered
  xspan <- mu_ext$xmax - mu_ext$xmin
  yspan <- mu_ext$ymax - mu_ext$ymin
  ncol  <- ceiling(xspan / desired_res_m)
  nrow  <- ceiling(yspan / desired_res_m)
  
  # Build grid anchored to Mukey xmin/ymin; expand if needed
  r <- terra::rast(xmin = mu_ext$xmin,
                   xmax = mu_ext$xmin + ncol * desired_res_m,
                   ymin = mu_ext$ymin,
                   ymax = mu_ext$ymin + nrow * desired_res_m,
                   ncols = ncol, nrows = nrow, crs = mu_crs)
  terra::values(r) <- NA
  r
}

# ---------------------------
# 1) AOI + buffer
# ---------------------------
# Read AOI; repair invalid geometries if needed; buffer in meters (assumes/projected CRS later).
aoi_sf <- sf::st_read(aoi_path, quiet = TRUE)
if (!all(sf::st_is_valid(aoi_sf))) {
  message("Repairing invalid AOI geometries...")
  aoi_sf <- sf::st_make_valid(aoi_sf)
}
aoi_buf <- sf::st_buffer(aoi_sf, buffer_m)

# ---------------------------
# 2) Mukey tiles ➜ UNION grid ➜ align (near) ➜ merge
# ---------------------------
# STRATEGY:
#   - To avoid huge WCS requests, tile the AOI buffer into a grid (tile_n).
#   - Request Mukey tiles per tile chunk at 30 m (categorical).
#   - Build a UNION raster extent using the first tile as the template, extending to cover others.
#   - Align each tile to the UNION with nearest neighbor (category-safe), then merge.
# NOTE:
#   - The UNION raster becomes the *anchor* for CRS, origin, and extent used later for final grid creation.
grid_mukey   <- sf::st_make_grid(aoi_buf, n = tile_n)
tiles_mukey  <- sf::st_intersection(aoi_buf, grid_mukey)
tiles_mukey  <- sf::st_collection_extract(tiles_mukey, "POLYGON")
aoi_list_mukey <- split(sf::st_as_sf(tiles_mukey), seq_len(nrow(tiles_mukey)))

mukey_tiles <- lapply(aoi_list_mukey, function(x) soilDB::mukey.wcs(aoi = x, db = "gNATSGO", res = 30))
stopifnot(length(mukey_tiles) > 0)

# Build UNION Mukey extent (keep CRS/res/origin from first tile)
mukey_union <- mukey_tiles[[1]]
if (length(mukey_tiles) > 1) {
  for (i in 2:length(mukey_tiles)) mukey_union <- terra::extend(mukey_union, mukey_tiles[[i]])
}

# Align each Mukey tile to UNION grid (nearest neighbor for categories)
mukey_aligned <- lapply(mukey_tiles, align_mukey_cat, target = mukey_union)

# Merge Mukey tiles (adjacent stitch) — EXPLICITLY use terra::merge
# NOTE:
#   - merge performs seamless mosaic for adjacent tiles with identical CRS/resolution.
mukey_template <- mukey_aligned[[1]]
if (length(mukey_aligned) > 1) {
  for (i in 2:length(mukey_aligned)) {
    mukey_template <- terra::merge(mukey_template, mukey_aligned[[i]])
  }
}
# Alternatively: mukey_template <- Reduce(terra::merge, mukey_aligned)

# ---------------------------
# 3) Create the FINAL Mukey-based grid at user resolution
#     (same CRS + origin + AOI extent; exact alignment preserved)
# ---------------------------
# RESULT:
#   - final_grid is empty (NA) but defines geometry; subsequent DEM tiles will be snapped to it.
final_grid <- make_final_grid_from_mukey(mukey_union, final_res_m)

plot(final_grid)

# Grid sanity — checks ensure we haven't drifted from the Mukey anchor.
stopifnot(terra::crs(final_grid) == terra::crs(mukey_template))
message(sprintf("Final DEM grid resolution: %.3f m", terra::res(final_grid)[1]))
message("Final DEM grid extent: ", as.character(terra::ext(final_grid)))
message("Mukey UNION extent:    ", as.character(terra::ext(mukey_union)))

# ---------------------------
# 4) Tile AOI ➜ fetch 3DEP ➜ normalize ➜ SNAP each DEM tile to FINAL grid
# ---------------------------
# Fetch directly in Mukey CRS so we avoid an extra reprojection step before snapping.
setCRS <- terra::crs(mukey_union)   # fetch directly in Mukey CRS
aoi_prj <- sf::st_transform(aoi_buf, setCRS)

# --- parameters ---
# Tile overlap ensures that after resampling/snap to final_grid, we don't lose edge pixels within tiles.
tile_overlap_m <- max(3 * final_res_m, 180)  # 3 cells or at least 180 m

# Build tile grid and add overlap per tile
tile_grid <- sf::st_make_grid(aoi_prj, n = tile_n)
tile_sf   <- sf::st_intersection(aoi_prj, tile_grid)
tile_sf   <- sf::st_collection_extract(tile_sf, "POLYGON")

# ADD OVERLAP
tile_sf_ol <- sf::st_buffer(tile_sf, tile_overlap_m)

# Keep overlap within the buffered AOI to avoid fetching far outside
tile_sf_ol <- sf::st_intersection(tile_sf_ol, aoi_prj)

# Convert to list for per-tile iterated fetching
tile_list <- split(sf::st_as_sf(tile_sf_ol), seq_len(nrow(tile_sf_ol)))

# Robust tile fetcher with retries
# NOTES:
#   - elevatr::get_elev_raster parameters:
#       - locations: sf polygon for clipping
#       - z: zoom level (native tile resolution)
#       - clip="locations": clip server-side to polygon
#       - prj: target CRS for the returned raster (Mukey CRS)
#       - src="aws": use AWS-hosted 3DEP tiles
#       - override_size_check: allow larger-than-default requests
#       - tmp_dir/ncpu: performance and local scratch management
#       - neg_to_na=FALSE: preserve negative elevations if present
get_tile <- function(poly_sf, z, prj, tmp_dir, tries = 3, sleep_sec = 3, ncpu = ncpu_fetch) {
  for (i in seq_len(tries)) {
    res <- try(
      elevatr::get_elev_raster(
        locations = poly_sf,
        z = z,
        clip = "locations",
        prj = prj,       # Mukey CRS
        src = "aws",
        override_size_check = TRUE,
        tmp_dir = tmp_dir,
        ncpu = ncpu,
        verbose = TRUE,
        neg_to_na = FALSE
      ),
      silent = TRUE
    )
    if (!inherits(res, "try-error")) return(res)
    message(sprintf("Tile attempt %d/%d failed; retrying in %ds ...", i, tries, sleep_sec))
    Sys.sleep(sleep_sec)
  }
  stop("Tile download failed after retries.")
}

message("Requesting 3DEP DEM tiles (src='aws', z=", dem_z, ") ...")
dem_tiles_raw <- purrr::map(
  tile_list,
  ~ get_tile(.x, z = dem_z, prj = setCRS, tmp_dir = tmp_dir,
             tries = retry_tries, sleep_sec = retry_backoff_s, ncpu = ncpu_fetch)
)

# Normalize NoData ➜ snap to FINAL grid ➜ write per-tile GeoTIFFs
# IMPORTANT:
#   - align_dem_bilinear(xi, final_grid) does both projection/resampling and grid snapping.
#   - Writing with LZW compression + tiling improves read-time performance for downstream analyses.
for (i in seq_along(dem_tiles_raw)) {
  xi <- dem_tiles_raw[[i]]
  if (!inherits(xi, "SpatRaster")) xi <- terra::rast(xi)
  xi <- normalize_nodata(xi)
  xi <- align_dem_bilinear(xi, final_grid)   # SNAP to user resolution on Mukey extent/CRS
  
  tile_path <- file.path(tmp_dir, sprintf("dem_tile_final_%02d.tif", i))
  terra::writeRaster(
    xi, filename = tile_path, overwrite = TRUE,
    wopt = list(
      gdal     = c("COMPRESS=LZW","TILED=YES","BIGTIFF=IF_NEEDED","BLOCKXSIZE=256","BLOCKYSIZE=256","PREDICTOR=3"),
      datatype = "FLT4S",
      NAflag   = -9999
    )
  )
}

# ---------------------------
# 5) Stack aligned DEM tiles ➜ aggregate with NA-safe mean ➜ write final
# ---------------------------
# The per-tile outputs now share identical geometry (final_grid); stacking is safe and deterministic.
dem.list <- list.files(path = tmp_dir, pattern = "^dem_tile_final_\\d+\\.tif$", full.names = TRUE)
stopifnot(length(dem.list) > 0)

tiles_final <- lapply(dem.list, terra::rast)  # identical geometry
dem_stack   <- terra::rast(tiles_final)       # multi-layer stack

message("Assembling final DEM: ")
# Aggregate per-cell using NA-safe mean across layers
# RATIONALE: mean reduces seam artifacts across adjacent tiles; alternative strategies include median or
#            weighted mosaics if overlaps vary; here, uniform overlaps + bilinear resampling make mean reasonable.
dem_final   <- terra::app(dem_stack, fun = mean_na)

# Assert final DEM grid equals the requested Mukey-based grid
# These checks enforce exact match on geometry; any mismatch indicates resampling/CRS drift.
stopifnot(terra::compareGeom(dem_final, final_grid, stopOnError = TRUE))
stopifnot(terra::ext(dem_final) == terra::ext(final_grid))
stopifnot(all(terra::res(dem_final) == terra::res(final_grid)))
stopifnot(terra::crs(dem_final) == terra::crs(final_grid))

# Stats & write
# minmax computes cached range; terra::global prints summary stats (ignoring NAs).
terra::minmax(dem_final)
print(terra::global(dem_final, c("min","max","mean"), na.rm = TRUE))

# Final write with LZW compression; BigTIFF fallback enabled for large outputs.
terra::writeRaster(
  dem_final, filename = out_final_dem, overwrite = TRUE,
  wopt = list(
    gdal     = c("COMPRESS=LZW","TILED=YES","BIGTIFF=IF_NEEDED","PREDICTOR=3"),
    datatype = "FLT4S",
    NAflag   = -9999
  )
)
message("Mukey-snapped DEM written: ", out_final_dem)

# ---------------------------
# 6) Diagnostics — verify NAflag & sentinels per tile
# ---------------------------
# PURPOSE:
#   - Confirm that tile outputs have the expected NA flag (-9999).
#   - Check for presence of extreme sentinel values (>1e20 or < -1e20) after normalization.
#   - Provide per-tile summary (min/max) for quick sanity checks.
diag <- lapply(seq_along(dem.list), function(i) {
  r <- terra::rast(dem.list[i])
  data.frame(
    tile     = basename(dem.list[i]),
    naflag   = suppressWarnings(tryCatch(terra::NAflag(r), error = function(e) NA)),
    min      = terra::global(r, "min",  na.rm = TRUE)[1,1],
    max      = terra::global(r, "max",  na.rm = TRUE)[1,1],
    gt1e20   = terra::global(terra::ifel(r >  1e20, 1, 0), "sum", na.rm = TRUE)[1,1],
    lt1e20   = terra::global(terra::ifel(r < -1e20, 1, 0), "sum", na.rm = TRUE)[1,1]
  )
})
print(do.call(rbind, diag))

# Final assertion: output grid equals FINAL grid
# This re-check protects against any file I/O surprises that might have altered metadata during write.
stopifnot(terra::compareGeom(terra::rast(out_final_dem), final_grid, stopOnError = TRUE))

message("Complete")
plot(dem_final)
