# ==========================================================
# Script to gather USGS 3DEP Best Available Elevation Data for a given
# area, mosaic, clip, return and save as a .tif
#
# Workflow: Load Packages -> load Boundary -> get elevation data -> save file
# Author: Dave White (with integrated enhancements)
# Date: 2026-01-26
# ==========================================================

suppressPackageStartupMessages({
  # Load required packages; install if missing
  required.packages <- c("elevatr", "terra", "sf")
  new.packages <- required.packages[!(required.packages %in% installed.packages()[, "Package"])]
  if (length(new.packages)) install.packages(new.packages)
  lapply(required.packages, require, character.only = TRUE)
})



#temp directory to store file downloads
dir <- "D:/federal_lands/26_Spodic_Intensity/temp"

# bring in your area of interest
aoi_sf <- st_read("D:/federal_lands/26_Spodic_Intensity/SpodicIntensityDSM_Boundary/SpodosolHuc12Opt2.shp")
aoi_buf <- st_buffer(aoi_sf, 100) # slight buffer
aoi_df <- data.frame(x = runif(6,min=sf::st_bbox(aoi_buf)$xmin, 
                               max=sf::st_bbox(aoi_buf)$xmax),
                     y = runif(6,min=sf::st_bbox(aoi_buf)$ymin, 
                               max=sf::st_bbox(aoi_buf)$ymax))

# 2. Get 3DEP Data using get_elev_raster
# z=14 is roughly 10m/1/3 arc-second. Use higher z (e.g., 16) for 1m if available.
# src = "aws" is generally faster and reliable for best available.
dem <- get_elev_raster(locations = aoi_df, 
                       z = 14, 
                       prj = st_crs(aoi_sf),
                       src = "aws",
                       override_size_check = T,
                       #tmp_dir = dir,
                       ncpu=detectCores()-1)

# 3. Convert to terra object for faster processing
dem_terra <- rast(dem)

# 4. Clip to specific polygon (optional, if you have a shapefile)
# If you only want the bounding box, step 2 already did that.
clipped_dem <- mask(crop(dem_terra, aoi_sf), aoi_sf)

# 5. Visualize
plot(clipped_dem, main = "USGS 3DEP Best Available Data")

# 6. Save the mosaiced/clipped DEM
writeRaster(clipped_dem, "D:/federal_lands/26_Spodic_Intensity/3dep_data_clipped.tif", overwrite=TRUE)
