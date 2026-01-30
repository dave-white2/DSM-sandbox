

#!/usr/bin/env Rscript
# ==========================================================
# Spodic Intensity Covariate Preparation
# Stack -> Tiles -> VRT (preserve categories) + Overviews + RAT fallback
# Author: Dave White + Integrated enhancements
# Date: 2026-01-26
# ==========================================================

suppressPackageStartupMessages({
  required.packages <- c("terra", "sf", "gdalUtilities", "fs", "stringr", "dplyr")
  new.packages <- required.packages[!(required.packages %in% installed.packages()[, "Package"])]
  if (length(new.packages))
    install.packages(new.packages)
  lapply(required.packages, require, character.only = TRUE)
})

## ---- 0) Environment --------------------------------------------------------
# Preserve factor levels via GDAL PAM (.aux.xml)
Sys.setenv(GDAL_PAM_ENABLED = "YES")

# GDAL cache: prefer numeric MB for portability (e.g., 4096 MB). If your GDAL supports % you can keep it.
# Sys.setenv(GDAL_CACHEMAX = "4096")   # MB
Sys.setenv(GDAL_CACHEMAX = "75%")      # Percent (ok if your GDAL supports it)

# Terra runtime options: use fast SSD tempdir if available
terraOptions(tempdir = "~/spodic/tiles/.terra_tmp",
             memfrac = 0.8,
             progress = 1)
if (!dir.exists("~/spodic/tiles/.terra_tmp"))
  dir.create("~/spodic/tiles/.terra_tmp", recursive = TRUE)

# Optional: log GDAL version
gdalinfo_bin <- Sys.which("gdalinfo")
if (nzchar(gdalinfo_bin)) {
  ver <- try(system2(gdalinfo_bin, "--version", stdout = TRUE), silent = TRUE)
  if (!inherits(ver, "try-error"))
    message("GDAL: ", paste(ver, collapse = " "))
}

## ---- 1) Paths --------------------------------------------------------------
path_spec <- "~/spodic/cov/spec"       # spectral (continuous)
path_cont <- "~/spodic/cov/ter/cont"   # terrain (continuous)
path_them <- "~/spodic/cov/ter/them"   # terrain (categorical)

out_dir    <- "~/spodic/tiles"
vrt_path   <- file.path(out_dir, "cov.vrt")
levels_rds <- file.path(out_dir, "categorical_levels.rds")   # cache of RATs for fallback

stopifnot(dir.exists(path_spec),
          dir.exists(path_cont),
          dir.exists(path_them))
if (!dir.exists(out_dir))
  dir.create(out_dir, recursive = TRUE)

## ---- 2) Load & align rasters ----------------------------------------------
clean_names <- function(paths)
  tools::file_path_sans_ext(basename(paths))

spec_files <- list.files(path_spec, pattern = "\\.tif$", full.names = TRUE)
cont_files <- list.files(path_cont, pattern = "\\.tif$", full.names = TRUE)
them_files <- list.files(path_them, pattern = "\\.tif$", full.names = TRUE)
stopifnot(length(spec_files) > 0,
          length(cont_files) > 0,
          length(them_files) > 0)

# Spectral + continuous (align spectral to the continuous grid)
rSPEC <- rast(spec_files)
names(rSPEC) <- clean_names(spec_files)
rCONT <- rast(cont_files)
names(rCONT) <- clean_names(cont_files)

if (crs(rSPEC) == crs(rCONT)) {
  rSPEC_aligned <- resample(rSPEC, rCONT, method = "bilinear")
  message("Spectral: resampled to continuous grid (bilinear)")
} else {
  rSPEC_aligned <- project(rSPEC, rCONT, method = "bilinear")
  message("Spectral: reprojected to continuous grid (bilinear)")
}
rCONT <- c(rCONT, rSPEC_aligned)

# Thematic (categorical): ensure factor type, align via nearest neighbor
rTHEM <- rast(them_files)
names(rTHEM) <- clean_names(them_files)
rTHEM <- as.factor(rTHEM)
if (!compareGeom(rTHEM, rCONT, stopOnError = FALSE)) {
  rTHEM <- project(rTHEM, rCONT, method = "near")
  message("Thematic: aligned to continuous grid (near)")
}
message("Is factor (rTHEM)? -> ", paste(is.factor(rTHEM), collapse = ", "))

# Full stack for tiling
rStack <- c(rCONT, rTHEM)
n_cont <- nlyr(rCONT)
n_them <- nlyr(rTHEM)
message(sprintf(
  "Stack dims: %d rows x %d cols x %d bands",
  nrow(rStack),
  ncol(rStack),
  nlyr(rStack)
))

## ---- 3) Cache factor levels (RAT) for fallback ----------------------------
# Store RATs from the categorical stack; used to re-apply after VRT if needed
lev_list <- levels(rTHEM)
saveRDS(lev_list, levels_rds)
message("Saved categorical factor levels to: ", levels_rds)

# Helper for RAT re-apply (expects lev_list[[i]] data.frames with at least 'ID')
apply_factor_levels <- function(x_cat, infile) {
  if (!file.exists(infile)) {
    warning("Factor-level cache not found: ",
            infile,
            " — levels will not be re-applied.")
    return(x_cat)
  }
  lev_list <- readRDS(infile)
  stopifnot(length(lev_list) == nlyr(x_cat))
  for (i in seq_len(nlyr(x_cat))) {
    x_cat[[i]] <- as.factor(x_cat[[i]])
    # Optional: verify 'ID' column presence
    if (!("ID" %in% names(lev_list[[i]]))) {
      warning("Level table for band ",
              i,
              " lacks 'ID' column; attaching anyway.")
    }
    levels(x_cat, i) <- lev_list[[i]]
  }
  x_cat
}

## ---- 4) Automatic tile sizing ---------------------------------------------
ram_gb_per_core <- 8
bytes_per_pixel <- 4
n_layers <- nlyr(rStack)

estimate_tile_side <- function(ram_gb,
                               n_layers,
                               bytes_per_pixel = 4,
                               safety_frac = 0.75) {
  usable_bytes <- ram_gb * 1e9 * safety_frac
  bytes_per_stack_pixel <- bytes_per_pixel * n_layers
  max_pixels <- usable_bytes / bytes_per_stack_pixel
  tile_side <- floor(sqrt(max_pixels))
  max(512L, tile_side)
}

tile_side <- estimate_tile_side(ram_gb_per_core, n_layers)
nr <- nrow(rStack)
nc <- ncol(rStack)
tile_rows <- ceiling(nr / tile_side)
tile_cols <- ceiling(nc / tile_side)
tile_grid <- rast(
  ncol = tile_cols,
  nrow = tile_rows,
  extent = ext(rStack),
  crs = crs(rStack)
)

message(
  sprintf(
    "Auto tile size ~%d px per side -> grid %d rows x %d cols",
    tile_side,
    tile_rows,
    tile_cols
  )
)

## ---- 5) Create tiles directly from stack ----------------------------------
setwd(out_dir)
message("Creating tiles directly from in-memory stack...")
tiles_info <- makeTiles(
  rStack,
  tile_grid,
  overwrite = TRUE,
  extend = FALSE,
  gdal = c(
    "TFW=YES",
    # optional; not strictly needed for GeoTIFF
    "COMPRESS=LZW",
    "TILED=YES",
    "BLOCKXSIZE=512",
    "BLOCKYSIZE=512",
    "BIGTIFF=IF_NEEDED"
  )      # IF_NEEDED is generally fine; YES also OK
)

# Gather a list of names from the tile GeoTIFFs
tl.names <- list.files(out_dir, pattern = "\\.tif$", full.names = TRUE)
stopifnot(length(tl.names) > 0)
message("Wrote ", length(tl.names), " tile(s).")

## ---- 6) Helper: robust gdaladdo wrapper via system2 -----------------------
# NOTE: gdaladdo requires '--config KEY VALUE' pairs, not 'KEY=VALUE'.
gdal_addo <- function(srcfile,
                      levels = c(2, 4, 8, 16),
                      resampler = c(
                        "nearest",
                        "average",
                        "rms",
                        "gauss",
                        "cubic",
                        "cubicspline",
                        "lanczos",
                        "average_mp",
                        "average_magphase",
                        "mode"
                      ),
                      ro = TRUE,
                      clean = FALSE,
                      config_options = c(
                        "COMPRESS_OVERVIEW=LZW",
                        "BIGTIFF_OVERVIEW=YES",
                        # force BigTIFF to avoid 4GB limit (safe for VRT .ovr)
                        "GDAL_TIFF_OVR_BLOCKSIZE=256",
                        "GDAL_NUM_THREADS=ALL_CPUS"
                      )) {
  stopifnot(file.exists(srcfile))
  resampler <- match.arg(resampler)
  
  gdal_bin <- Sys.which("gdaladdo")
  if (!nzchar(gdal_bin)) {
    warning("GDAL 'gdaladdo' not found on PATH. Skipping overviews for: ",
            srcfile)
    return(invisible(
      list(
        status = 1L,
        message = "gdaladdo not found",
        file = srcfile
      )
    ))
  }
  
  # Optionally clean existing overviews
  if (clean) {
    try(invisible(system2(
      gdal_bin,
      c("-clean", srcfile),
      stdout = TRUE,
      stderr = TRUE
    )), silent = TRUE)
  }
  
  # Expand config options: "KEY=VALUE" -> c("--config", "KEY", "VALUE")
  cfg_args <- character(0)
  if (length(config_options)) {
    kv <- strsplit(config_options, "=", fixed = TRUE)
    for (pair in kv) {
      if (length(pair) == 2) {
        cfg_args <- c(cfg_args, "--config", pair[[1]], pair[[2]])
      } else {
        warning("Skipping malformed config option: ",
                paste(pair, collapse = "="))
      }
    }
  }
  
  # Compose args: config, resampler, external overviews (-ro), file, levels
  args <- c(cfg_args,
            "-r",
            resampler,
            if (ro)
              "-ro"
            else
              NULL,
            srcfile,
            as.character(levels))
  
  out <- tryCatch({
    res <- system2(gdal_bin, args, stdout = TRUE, stderr = TRUE)
    status <- attr(res, "status")
    if (is.null(status))
      status <- 0L
    list(status = status, text = res)
  }, error = function(e)
    list(status = 1L, text = e$message))
  
  if (out$status == 0L) {
    message("     ✓ Overviews built for: ", basename(srcfile))
    invisible(list(
      status = 0L,
      message = "OK",
      file = srcfile
    ))
  } else {
    msg <- paste(out$text, collapse = "\n")
    message("     ✗ Overviews failed for ", basename(srcfile), ":\n", msg)
    invisible(list(
      status = 1L,
      message = msg,
      file = srcfile
    ))
  }
}

## ---- 7) (Optional) Build overviews per tile --------------------------------
# Building per-source overviews scales better than one huge VRT .ovr
message("Building GDAL overviews for all tiles (per-source)...")
resampler <- "nearest"  # consider "average" for continuous-only tiles
ovr_levels <- c(2, 4, 8, 16)

for (f in tl.names) {
  message("  ▶ ", basename(f))
  # Clean existing overviews and rebuild
  gdal_addo(
    srcfile = f,
    levels = ovr_levels,
    resampler = resampler,
    ro = TRUE,
    clean = TRUE,
    config_options = c(
      "COMPRESS_OVERVIEW=LZW",
      "BIGTIFF_OVERVIEW=IF_NEEDED",
      # per-source is typically fine with IF_NEEDED
      "GDAL_TIFF_OVR_BLOCKSIZE=256",
      "GDAL_NUM_THREADS=ALL_CPUS"
    )
  )
}

## ---- 8) Build VRT and (optionally) VRT overviews ---------------------------
message("Building VRT from tiles...")

build_vrt_with_overviews <- function(tl.names,
                                     vrt_path,
                                     resampler = "nearest",
                                     levels = c(2, 4, 8, 16),
                                     compress_overview = "LZW",
                                     ovr_blocksize = 256,
                                     add_vrt_overviews = TRUE) {
  stopifnot(is.character(tl.names), length(tl.names) > 0)
  stopifnot(is.character(vrt_path), nzchar(vrt_path))
  
  terra::vrt(
    x = tl.names,
    filename = vrt_path,
    overwrite = TRUE,
    set_names = TRUE
  )
  message("VRT built: ", vrt_path)
  
  if (!file.exists(vrt_path)) {
    msg <- sprintf("VRT path does not exist after build: %s", vrt_path)
    warning(msg)
    return(invisible(list(
      status = 1L,
      message = msg,
      vrt = vrt_path
    )))
  }
  
  if (!add_vrt_overviews) {
    return(invisible(
      list(
        status = 0L,
        message = "VRT built; overviews skipped",
        vrt = vrt_path
      )
    ))
  }
  
  # External overviews on VRT (BigTIFF forced to avoid 4GB limit)
  res <- gdal_addo(
    srcfile = vrt_path,
    levels = levels,
    resampler = resampler,
    ro = TRUE,
    clean = TRUE,
    config_options = c(
      paste0("COMPRESS_OVERVIEW=", compress_overview),
      "BIGTIFF_OVERVIEW=YES",
      # force BigTIFF on VRT .ovr
      paste0("GDAL_TIFF_OVR_BLOCKSIZE=", ovr_blocksize),
      "GDAL_NUM_THREADS=ALL_CPUS"
    )
  )
  if (res$status == 0L) {
    message("     ✓ VRT overviews built")
  }
  res
}

# Build VRT and choose whether to add overviews (TRUE = add; FALSE = rely on per-tile overviews)
res_vrt <- build_vrt_with_overviews(
  tl.names  = tl.names,
  vrt_path  = vrt_path,
  resampler = "nearest",
  # or "average" if prioritizing continuous bands
  levels    = c(2, 4, 8, 16),
  compress_overview = "LZW",
  add_vrt_overviews = F#TRUE       # set FALSE if you prefer only per-tile overviews
)

if (res_vrt$status != 0L) {
  warning("VRT/overview step reported issues: ", res_vrt$message)
}

## ---- 9) Validate categorical preservation on a sample tile ----------------
tile1 <- rast(tl.names[1])
message("Sample tile factor flags: ", paste(is.factor(tile1)[1:min(12, nlyr(tile1))], collapse = ", "))

## ---- 10) Load VRT & fallback RAT re-apply if needed -----------------------
rStack_vrt <- rast(vrt_path)
message("Bands in VRT: ", nlyr(rStack_vrt))

# Check factor flags for THEM portion
factor_flags <- is.factor(rStack_vrt)
message("Factor flags for THEM bands: ",
        paste(factor_flags[(n_cont + 1):(n_cont + n_them)], collapse = ", "))

# If any THEM bands lost factor status, re-apply cached RAT
if (!all(factor_flags[(n_cont + 1):(n_cont + n_them)])) {
  message("Re-applying categorical factor levels to VRT from cache...")
  them_vrt <- rStack_vrt[[(n_cont + 1):(n_cont + n_them)]]
  them_vrt <- apply_factor_levels(them_vrt, levels_rds)
  rStack_vrt <- c(rStack_vrt[[1:n_cont]], them_vrt)  # replace THEM portion in-memory
  message("After re-apply, factor flags (THEM): ",
          paste(is.factor(rStack_vrt)[(n_cont + 1):(n_cont + n_them)], collapse = ", "))
}

message("Done.")
