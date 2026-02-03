# ==========================================================
# Script: USGS 3DEP Best Available Elevation for AOI
# Workflow: Load Packages -> Load Boundary -> Get DEM -> Smooth (optional)
#           -> Align to Mukey template -> Clip/Mask -> Save GeoTIFF
# Author: Dave White
# Date: 2026-01-26 (revised)
# ==========================================================

suppressPackageStartupMessages({
  required.packages <- c("elevatr", "terra", "sf", "parallel", "soilDB")
  new.packages <- required.packages[!(required.packages %in% installed.packages()[, "Package"])]
  if (length(new.packages)) install.packages(new.packages)
  lapply(required.packages, require, character.only = TRUE)
})

# ---------------------------
# Paths & parameters
# ---------------------------
tmp_dir <- "D:/federal_lands/26_Spodic_Intensity/temp"
aoi_path <- "D:/federal_lands/26_Spodic_Intensity/SpodicIntensityDSM_Boundary/SpodosolHuc12Opt2.shp"
out_tif  <- "D:/federal_lands/26_Spodic_Intensity/3dep_data_clipped.tif"

# Buffer distance (meters) around AOI
buffer_m <- 500

# Smoothing radius (meters) prior to resampling (set to 0 to skip)
smooth_radius_m <- 20

# Target DEM zoom level; z=14 ~10 m; z=16 ~1 m (if available)
dem_z <- 14

# ---------------------------
# 1. Load AOI and buffer
# ---------------------------
aoi_sf <- st_read(aoi_path, quiet = TRUE)
if (!all(st_is_valid(aoi_sf))) {
  message("Repairing invalid AOI geometries...")
  aoi_sf <- st_make_valid(aoi_sf)
}
aoi_buf <- st_buffer(aoi_sf, buffer_m)

# ---------------------------
# 2. Create Mukey template using soilDB (gNATSGO, ~30m)
#    This gives us the target CRS and grid for final alignment.
# ---------------------------
# Split AOI into 2x2 tiles to avoid large requests for WCS
grid <- st_make_grid(aoi_buf, n = c(2, 2))
plot(grid)
tiles <- st_intersection(aoi_buf, grid)
plot(tiles$geometry)
tiles <- st_combine(tiles)
plot(tiles)
single_polygons_sf <- st_cast(tiles, "POLYGON", warn = FALSE)
plot(single_polygons_sf)
single_polygons_sf <- st_as_sf(single_polygons_sf)
aoi_list <- as.list(split(single_polygons_sf, seq_len(nrow(single_polygons_sf))))
plot(aoi_list$`1`)

raster_list <- lapply(aoi_list, function(x) {
  # mukey.wcs returns SpatRaster when terra is loaded
  soilDB::mukey.wcs(aoi = x, db = "gNATSGO", res = 30)
})

# Merge Mukey tiles -> Template
merged_mukey <- terra::mosaic(terra::sprc(raster_list), fun="mean")
plot(merged_mukey, main = "Mukey template (gNATSGO)", col = hcl.colors(40, "viridis"), axes = FALSE)

# Template CRS (use this as the DEM target CRS)
setCRS <- terra::crs(merged_mukey)
stopifnot(!is.na(setCRS))

# Reproject AOI to template CRS
aoi_prj <- st_transform(aoi_buf, setCRS)

# ---------------------------
# 3. Fetch 3DEP Best Available DEM via elevatr (AWS)
#    We pass the AOI in the Mukey CRS to avoid reprojection later.
# ---------------------------
message("Requesting 3DEP DEM (z=", dem_z, ", src='aws') ...")
dem_raster <- {
  elevatr::get_elev_raster(
    locations = aoi_prj,
    z = dem_z,
    clip = "locations",
    prj = setCRS,
    src = "aws",
    override_size_check = TRUE,
    tmp_dir = tmp_dir,
    ncpu = parallel::detectCores() - 1,
    verbose = T,
    neg_to_na = T
  )
}

plot(dem_raster)

dem_terra <- terra::rast(dem_raster)

plot(dem_terra)

# ---------------------------
# 7. Write GeoTIFF with sensible GDAL options
# ---------------------------
# Keep float precision for smoothed/resampled elevation
writeRaster(
  dem_terra,
  filename = out_tif,
  overwrite = TRUE,
  filetype = "GTiff",
  gdal = c("COMPRESS=LZW", "PREDICTOR=3", "TILED=YES", "BIGTIFF=IF_SAFER")
)
message("Saved: ", out_tif)

dem_terra <- terra::rast("D:/federal_lands/26_Spodic_Intensity/3dep_data_clipped.tif")

# check resolution of dem
dem_terra

# Resample DEM to 10m grid (bilinear for continuous elevation)
dem_terra <- terra::resample(dem_terra, 3, method = "bilinear", threads=T)
# check resolution of dem
dem_terra
plot(dem_terra)

# ---------------------------
# 5. Optional smoothing before alignment/resampling
# ---------------------------
dem_for_smooth <- dem_terra
if (smooth_radius_m > 0) {
  message("Smoothing DEM with circular window, radius=", smooth_radius_m, " m")
  # Build focal window in map units using terra::focalMat
  fw <- terra::focalMat(dem_for_smooth, d = smooth_radius_m, type = "circle")
  # Use a binary mask to apply mean over the neighborhood
  fw_binary <- fw
  fw_binary[fw_binary > 0] <- 1
  dem_for_smooth <- terra::focal(dem_for_smooth, w = fw_binary, fun = "mean")
}
plot(dem_for_smooth)


# ---------------------------
# 6. Align DEM to Mukey template (CRS & resolution), then clip/mask to AOI
# ---------------------------
# Project (if needed) — dem was already requested in setCRS, but keep this safe:
if (terra::crs(dem_for_smooth) != terra::crs(merged_mukey)) {
  dem_proj <- terra::project(dem_for_smooth, merged_mukey)
} else {
  dem_proj <- dem_for_smooth
}

# Resample DEM to Mukey grid (bilinear for continuous elevation)
dem_resamp <- terra::resample(dem_proj, merged_mukey, method = "bilinear", threads=T)

# Clip to AOI extent and mask to AOI polygon
# Use vector AOI in template CRS
clipped_dem <- terra::mask(terra::crop(dem_resamp, aoi_prj), aoi_prj)


# Stats & sanity visualization
plot(clipped_dem,
     main = "USGS 3DEP (aligned to Mukey, clipped to AOI)",
     col = hcl.colors(20, "Terrain"), axes = FALSE)

# ---------------------------
# 7. Write GeoTIFF with sensible GDAL options
# ---------------------------
# Keep float precision for smoothed/resampled elevation
writeRaster(
  clipped_dem,
  filename = out_tif,
  overwrite = TRUE,
  filetype = "GTiff",
  datatype = "FLT4S",
  gdal = c("COMPRESS=LZW", "PREDICTOR=3", "TILED=YES", "BIGTIFF=IF_SAFER")
)

message("Saved: ", out_tif)

# ---------------------------
# 8. Diagnostics (optional)
# ---------------------------
message("DEM datatype: ", terra::datatype(dem_terra))
message("DEM NAflag: ", paste0(terra::NAflag(dem_terra)))
print(terra::global(clipped_dem, fun = c("min", "mean", "max"), na.rm = TRUE))

