# ==========================================================
# Spodic Intensity Covariate Preparation
# Workflow: Stack → Data Reduction → Tiles → VRT + Overviews
# Author: Dave White (with integrated enhancements)
# Date: 2026-01-26
# ==========================================================

suppressPackageStartupMessages({
  # Load required packages; install if missing
  required.packages <- c("terra", "sf", "gdalUtilities", "fs", "stringr", "dplyr",
                         "parallel", "doParallel", "caret", "corrplot", "gdalUtilities", "xml2")
  new.packages <- required.packages[!(required.packages %in% installed.packages()[, "Package"])]
  if (length(new.packages)) install.packages(new.packages)
  lapply(required.packages, require, character.only = TRUE)
})

## ---- 0) Environment --------------------------------------------------------
# Enable GDAL PAM (.aux.xml) for factor level preservation
Sys.setenv(GDAL_PAM_ENABLED = "YES")

# Configure GDAL cache size (use MB or % depending on GDAL support)
# Sys.setenv(GDAL_CACHEMAX = "4096")   # Example: 4096 MB
Sys.setenv(GDAL_CACHEMAX = "75%")      # Example: 75% of available memory

# Set Terra runtime options (prefer SSD for temp directory)
terraOptions(tempdir = "~/spodic/tiles/.terra_tmp",
             memfrac = 0.8,
             progress = 1)

## ---- 1) Paths --------------------------------------------------------------
# Define input and output paths
path_spec <- "~/spodic/cov/spec"       # Spectral covariates (continuous)
path_ter  <- "~/spodic/cov/ter/cont"   # Terrain covariates (continuous)
out_dir   <- "~/spodic/tiles"          # Output directory for tiles
vrt_path  <- file.path(out_dir, "cov.vrt")  # VRT output path

## ---- 2) Load & align rasters ----------------------------------------------
# Helper: remove file extensions for clean layer names
clean_names <- function(paths) tools::file_path_sans_ext(basename(paths))

# List raster files
spec_files <- list.files(path_spec, pattern = "\\.tif$", full.names = TRUE)
ter_files  <- list.files(path_ter,  pattern = "\\.tif$", full.names = TRUE)

# Load rasters and assign names
rSPEC <- rast(spec_files); names(rSPEC) <- clean_names(spec_files)
rTER  <- rast(ter_files);  names(rTER)  <- clean_names(ter_files)

# ------------------------------
# Robust alignment + stacking
# ------------------------------
library(terra)

align_and_stack <- function(rSPEC, rTER, prefer = c("smaller", "rSPEC", "rTER")) {
  prefer <- match.arg(prefer)
  
  # --- 0) Sanity checks ---
  stopifnot(inherits(rSPEC, "SpatRaster"), inherits(rTER, "SpatRaster"))
  if (nlyr(rSPEC) == 0 || nlyr(rTER) == 0) stop("Both inputs must have at least one layer.")
  
  # --- 1) Choose a provisional template (by preference) ---
  template_pre <- switch(prefer,
                         "smaller" = if (ncell(rSPEC) <= ncell(rTER)) rSPEC else rTER,
                         "rSPEC"   = rSPEC,
                         "rTER"    = rTER
  )
  
  # --- 2) Harmonize CRS (only if needed) ---
  resample_method <- function(x) if (any(is.factor(x))) "near" else "bilinear"
  if (!identical(crs(rSPEC), crs(template_pre))) {
    rSPEC <- project(rSPEC, template_pre, method = resample_method(rSPEC))
  }
  if (!identical(crs(rTER), crs(template_pre))) {
    rTER  <- project(rTER,  template_pre, method = resample_method(rTER))
  }
  
  # --- 3) Crop both to the overlapping extent (intersection) ---
  e <- ext(
    max(xmin(rSPEC), xmin(rTER)),
    min(xmax(rSPEC), xmax(rTER)),
    max(ymin(rSPEC), ymin(rTER)),
    min(ymax(rSPEC), ymax(rTER))
  )
  # Explicit width/height checks (SpatExtent-safe)
  if ((xmax(e) - xmin(e)) <= 0 || (ymax(e) - ymin(e)) <= 0) {
    stop("No overlapping area after CRS alignment; cannot stack.")
  }
  rSPEC <- crop(rSPEC, e, snap = "out")
  rTER  <- crop(rTER,  e, snap = "out")
  
  # --- 4) Build a template grid on template_pre's grid, cropped to 'e' ---
  # This preserves resolution & origin (grid), but limits extent to the intersection.
  template <- crop(template_pre, e, snap = "out")
  
  # --- 5) Resample both stacks EXACTLY to the template grid ---
  resample_to_template <- function(x, template) {
    idx_fac <- which(is.factor(x))  # logical per layer; which() gets indices
    if (length(idx_fac) == 0) {
      out <- resample(x, template, method = "bilinear")
      names(out) <- names(x)
      out
    } else if (length(idx_fac) == nlyr(x)) {
      out <- resample(x, template, method = "near")
      names(out) <- names(x)
      out
    } else {
      # Mixed: resample numeric and factor separately, then reassemble in original order
      idx_num <- setdiff(seq_len(nlyr(x)), idx_fac)
      x_num <- resample(x[[idx_num]], template, method = "bilinear")
      names(x_num) <- names(x)[idx_num]
      x_fac <- resample(x[[idx_fac]], template, method = "near")
      names(x_fac) <- names(x)[idx_fac]
      out <- c(x_num, x_fac)
      out <- out[[match(names(x), names(out))]]  # restore original layer order
      out
    }
  }
  
  rSPEC2 <- resample_to_template(rSPEC, template)
  rTER2  <- resample_to_template(rTER,  template)
  
  # --- 6) Final validation ---
  # Terra 1.8.93: compareGeom does not accept named flags; keep it simple and
  # fall back to a manual check if compareGeom errors.
  ok <- tryCatch({
    terra::compareGeom(rSPEC2, rTER2)  # returns TRUE or throws
  }, error = function(e) NA)
  
  if (is.na(ok)) {
    # Manual geometry check with tolerance
    tol <- 1e-9
    same_crs    <- identical(crs(rSPEC2), crs(rTER2))
    same_res    <- all(abs(res(rSPEC2)    - res(rTER2))    < tol)
    same_origin <- all(abs(origin(rSPEC2) - origin(rTER2)) < tol)
    same_ext    <- all(abs(c(xmin(rSPEC2) - xmin(rTER2),
                             xmax(rSPEC2) - xmax(rTER2),
                             ymin(rSPEC2) - ymin(rTER2),
                             ymax(rSPEC2) - ymax(rTER2))) < tol)
    if (!all(c(same_crs, same_res, same_origin, same_ext))) {
      bad <- c("crs","res","origin","extent")[!c(same_crs, same_res, same_origin, same_ext)]
      stop("Geometry mismatch remains after alignment: ", paste(bad, collapse = ", "))
    }
  } else if (!isTRUE(ok)) {
    stop("Geometry mismatch remains after alignment.")
  }
  
  # --- 7) Stack ---
  stacked <- c(rSPEC2, rTER2)
  return(stacked)
}

# Combine into a single stack (prefer the smaller grid after intersection)
rStack <- align_and_stack(rSPEC, rTER, prefer = "smaller")

message(sprintf("Stack dimensions: %d rows x %d cols x %d bands",
                nrow(rStack), ncol(rStack), nlyr(rStack)))

# ------------------------------
# Helper: Standardize layer names
# ------------------------------
standardize_names <- function(x, max_len = 6) {
  y <- tolower(x)
  y <- gsub("_", "", y)
  y <- gsub("2021-", "", y)
  y <- gsub("-22", "", y)
  y <- gsub("-", "", y)
  y <- abbreviate(y, minlength = max_len, method = "both.sides")
  return(y)
}

# Apply standardized names to rStack
orig_names_all <- names(rStack)
std_names_all  <- standardize_names(orig_names_all, max_len = 6)

# Keep a mapping for reference
name_map <- data.frame(original = orig_names_all,
                       standardized = std_names_all,
                       stringsAsFactors = FALSE)
write.csv(name_map, file = "~/spodic/tiles/covNames.csv")
names(rStack) <- std_names_all

## ---- 3) Data Reduction -----------------------------------------------------
# Sample raster stack to a data frame
df <- spatSample(rStack, method = "regular", na.rm = TRUE, as.df = TRUE, size = 1000000)

# Verify column names match raster names
names(df)
names(rStack)

# ---- Near-zero variance filtering ------------------------------------------
# Remove covariates with near-zero variance (parallelized)
cl <- makeCluster(detectCores() - 1) # Use all but one core
registerDoParallel(cl)
zeroVar <- nearZeroVar(df, foreach = TRUE, allowParallel = TRUE)
stopCluster(cl)

gc() # Clean up memory

# Inspect zeroVar result (integer(0) means no removal needed)
head(zeroVar)

# Remove near-zero variance covariates
df2 <- if(length(zeroVar) > 0) df[, -zeroVar] else df

# Compare original vs reduced covariate count
length(df)
length(df2)

# Subset raster stack to reduced covariates
rStack <- subset(rStack, names(df2))

# Resample after reduction
df3 <- spatSample(rStack, method = "regular", na.rm = TRUE, as.df = TRUE, size = 1000000)

# ---- Correlation filtering -------------------------------------------------
corMat <- cor(df3)
corrplot(corMat, method = "circle")
highCorr <- findCorrelation(corMat, cutoff = 0.98)

df_reduced <- if (length(highCorr) > 0) df3[, -highCorr, drop = FALSE] else df3

# Subset raster stack again
rStack <- subset(rStack, names(df_reduced))

## ---- Recursive Feature Elimination (RFE) -------------------------------
# Load training dataset
train <- read_sf("~/spodic/data/train.shp")

# Keep only Spodic Index (SI)
train <- train["SI"]

# Extract covariate values for training points
pts.sv <- terra::extract(rStack, train, xy = FALSE, bind = TRUE, na.rm = TRUE, raw = FALSE)

# Convert to sf and drop geometry
comp <- st_as_sf(pts.sv)
comp2 <- na.omit(comp)
comp2 <- sf::st_drop_geometry(comp2)

# Prepare subsets for RFE
subsets <- seq(0, length(names(rStack)), 5)

# Cross-validation settings
number <- 10
repeats <- 5

# Set seeds for reproducibility
set.seed(12)
seeds <- vector(mode = "list", length = (number * repeats + 1))
for(i in 1:(number * repeats + 1)) seeds[[i]] <- sample.int(1000, number * repeats + 2)
seeds[[(number * repeats + 1)]] <- sample.int(1000, 1)

# Parallel setup for RFE
cl <- makeCluster(60) # Adjust cores; leave some for OS
registerDoParallel(cl)

# Configure RFE control
ctrl.RFE <- rfeControl(functions = rfFuncs,
                       method = "repeatedcv",
                       number = number,
                       repeats = repeats,
                       seeds = seeds,
                       verbose = FALSE)

# Run RFE
set.seed(9)
rfe <- rfe(x = comp2[, -c(1)],
           y = comp2$SI,
           sizes = subsets,
           method = "rf",
           rfeControl = ctrl.RFE,
           allowParallel = TRUE,
           metric = "RMSE")

stopCluster(cl)
gc()

# Inspect RFE results
rfe
plot(rfe)
predictors(rfe)

# Subset raster stack by selected predictors
rStack <- subset(rStack, predictors(rfe))

## ---- 4) Automatic Tile Sizing ---------------------------------------------
# Detect system resources
num_cores <- parallel::detectCores()
cat(sprintf("Number of cores: %d\n", num_cores))

# Get total RAM in GB
mem_info <- system("grep MemTotal /proc/meminfo", intern = TRUE)
mem_kb   <- as.numeric(strsplit(mem_info, "\\s+")[[1]][2])
mem_gb   <- mem_kb / 1024 / 1024
cat(sprintf("Total RAM in GB: %.2f\n", mem_gb))

# Estimate tile size based on RAM and stack size
ram_per_core_gb <- mem_gb / num_cores
ram_gb_per_core <- round(ram_per_core_gb)
bytes <- 8L  # Bytes per pixel (4 for float, 8 for double)

estimate_tile_side <- function(ram_gb, n_layers, bytes_per_pixel = bytes, safety_frac = 0.75) {
  usable_bytes <- ram_gb * 1e9 * safety_frac
  bytes_per_stack_pixel <- as.numeric(bytes_per_pixel) * n_layers
  max_pixels <- usable_bytes / bytes_per_stack_pixel
  tile_side <- floor(sqrt(max_pixels))
  max(512L, tile_side)
}

n_layers <- nlyr(rStack)
tile_side <- estimate_tile_side(ram_gb = ram_gb_per_core,
                                n_layers = n_layers,
                                bytes_per_pixel = bytes,
                                safety_frac = 0.75)
cat(sprintf("Estimated tile side (pixels): %d\n", tile_side))

# Create tile grid
nr <- nrow(rStack)
nc <- ncol(rStack)
tile_rows <- ceiling(nr / tile_side)
tile_cols <- ceiling(nc / tile_side)
tile_grid <- rast(ncol = tile_cols,
                  nrow = tile_rows,
                  extent = ext(rStack),
                  crs = crs(rStack))

message(sprintf("Auto tile size ~%d px per side → grid %d rows x %d cols",
                tile_side, tile_rows, tile_cols))

## ---- 5) Create Tiles -------------------------------------------------------
setwd(out_dir)
message("Creating tiles from in-memory stack...")
tiles_info <- makeTiles(rStack, tile_grid, overwrite = TRUE, extend = FALSE,
                        gdal = c("TFW=YES", "COMPRESS=LZW", "TILED=YES",
                                 "BLOCKXSIZE=512", "BLOCKYSIZE=512", "BIGTIFF=YES"))

# Verify tiles
tl.names <- list.files(out_dir, pattern = "\\.tif$", full.names = TRUE)
stopifnot(length(tl.names) > 0)
message("Wrote ", length(tl.names), " tile(s).")

# Inspect one tile's band descriptions as Terra sees them
one_tile <- tl.names[1]
tile_r <- terra::rast(one_tile)
print(names(tile_r))



## ---- 6) Build Overviews for Tiles ----------------------------------------
resampler  <- "average"
ovr_levels <- c(2, 4, 8, 16)

# Named character vector: names are the keys, values are strings
config_opts <- c(
  COMPRESS_OVERVIEW      = "LZW",
  BIGTIFF_OVERVIEW       = "IF_NEEDED",
  GDAL_TIFF_OVR_BLOCKSIZE= "256",      # quote numeric values as strings
  GDAL_NUM_THREADS       = "ALL_CPUS"
)

message("Building GDAL overviews for all tiles...")
for (f in tl.names) {
  message("  ▶ ", basename(f))
  gdal_addo(
    file           = f,
    overviews      = ovr_levels,
    method         = resampler,
    read_only      = TRUE,    # external .ovr for GeoTIFFs
    clean          = TRUE,    # True deletes.ovr files
    config_options = config_opts
  )
}
for (f in tl.names) {
  message("  ▶ ", basename(f))
  gdal_addo(
    file           = f,
    overviews      = ovr_levels,
    method         = resampler,
    read_only      = TRUE,    # external .ovr for GeoTIFFs
    clean          = F,       # F allows .ovr files to be created
    config_options = config_opts
  )
}

## ---- 7) Build VRT and Overviews -------------------------------------------
message("Building VRT from tiles...")
build_vrt_with_overviews <- function(tl.names, vrt_path,
                                     resampler = "average",
                                     levels = c(2, 4, 8, 16),
                                     compress_overview = "LZW",
                                     add_vrt_overviews = FALSE) {
  # Expand and normalize output path; ensure directory is writable
  vrt_path <- normalizePath(path.expand(vrt_path), winslash = "/", mustWork = FALSE)
  out_dir  <- dirname(vrt_path)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  stopifnot(file.access(out_dir, 2) == 0)
  
  # Absolute input paths
  tl.names <- normalizePath(tl.names, winslash = "/", mustWork = TRUE)
  
  # Build the mosaic VRT
  gdalbuildvrt(
    gdalfile   = tl.names,
    output.vrt = vrt_path,
    overwrite  = TRUE
  )  # gdalbuildvrt builds a VRT from a list of datasets. [3](https://community.pix4d.com/t/warning-w9020/6061)[4](https://www.rdocumentation.org/packages/sf/versions/1.0-23/topics/gdal_addo)
  
  # Optionally add overviews to the VRT (external .ovr)
  if (add_vrt_overviews) {
    sf::gdal_addo(
      file           = vrt_path,
      overviews      = levels,
      method         = resampler,
      read_only      = TRUE,  # external sidecar .ovr, equivalent to -ro in GDAL. [3](https://community.pix4d.com/t/warning-w9020/6061)
      clean          = FALSE,
      config_options = c(COMPRESS_OVERVIEW = compress_overview)  # named vector required. [2](https://www.rdocumentation.org/packages/gdalUtilities/versions/1.2.5/topics/gdalbuildvrt)
    )
  }
  return(vrt_path)
}


res_vrt <- build_vrt_with_overviews(
  tl.names = tl.names,
  vrt_path = file.path(path.expand("~/spodic/tiles"), "cov.vrt"),
  resampler = "average",
  levels = c(2, 4, 8, 16),
  compress_overview = "LZW",
  add_vrt_overviews = T
)

# 1) Use the names you want to persist
band_names <- names(rStack)        # or names(rStack_vrt) if you already subset
stopifnot(length(band_names) == terra::nlyr(rStack_vrt))

# 2) Read the VRT XML
doc <- read_xml(vrt_path)

# 3) Find all band nodes
bands <- xml_find_all(doc, "//VRTRasterBand")
stopifnot(length(bands) == length(band_names))

# 4) For each band, add (or replace) the <Description> tag
for (i in seq_along(bands)) {
  # Remove an existing Description if present
  desc_old <- xml_find_first(bands[[i]], "Description")
  if (!is.na(desc_old)) xml_remove(desc_old)
  
  # Add the desired description
  xml_add_child(bands[[i]], "Description", band_names[i])
}

# 5) Write back the VRT
write_xml(doc, vrt_path)


## ---- 9) Load VRT ---------------------------------------------------------
rStack_vrt <- rast(vrt_path)
message("Bands in VRT: ", nlyr(rStack_vrt))

names(rStack_vrt)


# Optional: quick grid diagnostics for debugging (comment/uncomment as needed)
# grid_diag <- function(x, tag = "") {
#   cat(sprintf(
#     "%s: crs=%s | res=(%.6f, %.6f) | origin=(%.6f, %.6f) | extent=[%.3f, %.3f, %.3f, %.3f] | n=(%d, %d, %d)\n",
#     tag, crs(x), res(x)[1], res(x)[2],
#     origin(x)[1], origin(x)[2],
#     xmin(x), xmax(x), ymin(x), ymax(x),
#     nrow(x), ncol(x), nlyr(x)
#   ))
# }
# grid_diag(rSPEC, "SPEC pre")
# grid_diag(rTER,  "TER  pre")
# grid_diag(rStack, "STACK")

