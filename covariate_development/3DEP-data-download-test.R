# ============================================================
# Get 3DEP-data
#
# Dave White - 2/4/26
# ============================================================
suppressPackageStartupMessages({
  pkgs <- c("elevatr","terra","sf")
  new  <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
  if (length(new)) install.packages(new)
  lapply(pkgs, require, character.only = TRUE)
})

# ---------------------------
# Paths & parameters (EDIT THESE)
# ---------------------------
tmp_dir         <- "D:/federal_lands/26_Spodic_Intensity/temp3"
aoi_path        <- "D:/federal_lands/26_Spodic_Intensity/test_poly.shp"
out_final_dem   <- "D:/federal_lands/26_Spodic_Intensity/dem_testPoly.tif"

dem_z           <- 12        # elevatr zoom (~5–10 m source tiles) — higher z => finer native tiles, larger downloads

ncpu            <- max(1, parallel::detectCores() - 2)  # parallel workers used by elevatr

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

# bring in boundary layer
poly_sf <- st_read(aoi_path)

# 
elev_dat <- elevatr::get_elev_raster(
  locations = poly_sf,
  z = dem_z,
  clip = "locations",
  src = "aws",
  override_size_check = TRUE,
  tmp_dir = tmp_dir,
  ncpu = ncpu,
  verbose = TRUE,
  neg_to_na = FALSE)


plot(elev_dat)


terra::writeRaster(
  elev_dat, filename = out_final_dem, overwrite = TRUE,
  wopt = list(
    gdal     = c("COMPRESS=LZW","TILED=YES","BIGTIFF=IF_NEEDED","BLOCKXSIZE=256","BLOCKYSIZE=256","PREDICTOR=3"),
    datatype = "FLT4S",
    NAflag   = -9999))
