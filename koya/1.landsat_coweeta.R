### Landsat 8 covariate development
# For the 2020 Florida Field Week

# 11/08/2019
# Dave White USDA-NRCS david.white@usda.gov
# Nicholas Klein-Baer USDA-USFS nicholas.klein-baer@usda.gov
# 

# The landsat data was dounloaded using the USFS implementations in google earth engine. The following is a link to those modules : https://earthengine.googlesource.com/users/USFS_GTAC/modules/+/master. To represent off periods of time Jan 1 to Feb 28th was used from 2015 to 2019. To represent wet periods of time July 1st to Sept 30th was used from 2015 to 2019. Time buffer was set to 2, weights was set to 1,2,3,2,1, median and SR was selected for each season. This algorithim calculates the median values for each date range. The final band combonations reflect those of Landsat 7.




## Load packages
required.packages <- c("raster", "sp", "rgdal", "RStoolbox", 
                       "snow", "snowfall","parallel", "itertools","doParallel")
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(required.packages, require, character.only=T)
rm(required.packages, new.packages)


# If using a laptop this section of code may be usefull for memory handling
## Increase actuve memory useable by raster package
#memory.limit(110000)
#rasterOptions(maxmemory = 5e+08, chunksize = 3e+07)

# Bring in our LS raster stacks and clip to our covariate development boundary

# load ls raster stacks

# set directory that contains the ls rasters 
# *** pay attention to the / below
setwd("C:/cov/ls8")


# read in raster layers to create raster stack
ls.off <- stack("leafoff.tif")


# visually inspect the raster
plot(ls.off$leafoff_1)



#### Extract individual bands from raster stacks.
# added a constant to each band. there were some issues in calculating the indices, where it was 
# returning no data cells. 


setwd("C:/cov/ls8/off")
lsoffb1 <- ls.off$leafoff_1 + 1
writeRaster(lsoffb1, filename = "lsoff1.tif", format = "GTiff", overwrite = TRUE)
lsoffb2 <- ls.off$leafoff_2 + 1
writeRaster(lsoffb2, filename = "lsoff2.tif", format = "GTiff", overwrite = TRUE)
lsoffb3 <- ls.off$leafoff_3 + 1
writeRaster(lsoffb3, filename = "lsoff3.tif", format = "GTiff", overwrite = TRUE)
lsoffb4 <- ls.off$leafoff_4 + 1
writeRaster(lsoffb4, filename = "lsoff4.tif", format = "GTiff", overwrite = TRUE)
lsoffb5 <- ls.off$leafoff_5 + 1
writeRaster(lsoffb5, filename = "lsoff5.tif", format = "GTiff", overwrite = TRUE)
lsoffb6 <- ls.off$leafoff_6 + 1
writeRaster(lsoffb6, filename = "lsoff6.tif", format = "GTiff", overwrite = TRUE)
lsoffb7 <- ls.off$leafoff_7 + 1
writeRaster(lsoffb7, filename = "lsoff7.tif", format = "GTiff", overwrite = TRUE)



######################################################################




## Development of Landsat off covariates

## set working directory for landsat 8 off data
setwd("C:/cov/ls8/off")

## load ls drt data as a raster stack
r.list=list.files(getwd(), pattern="tif$", full.names = FALSE)
ls.off <- stack(r.list)

## get individual bands
b1off <- ls.off$layer.1
b2off <- ls.off$layer.2
b3off <- ls.off$layer.3
b4off <- ls.off$layer.4
b5off <- ls.off$layer.5
b6off <- ls.off$layer.6
b7off <- ls.off$layer.7

## Normalized Difference index function
nd_fn <- function(bd1,bd2) {ind <- (bd1 - bd2)/(bd1 + bd2)*100
return(ind)
}

## set up cluster
beginCluster()
#### set workspace for leaf off covariates
setwd("C:/cov/ls8/off")


## ratio calcs for off
# note the  compression and datatypes are commented out, check the data ranges of each raster produced and select the appropriate datatype. Use ?raster::datatype to see the different choices.
s3t1off <- stack(b3off,b1off) 
clusterR(s3t1off, overlay, args=list(fun=nd_fn),progress='text',filename="nd3t1off.tif", options=c("COMPRESS=DEFLATE"),datatype='INT2S', overwrite=TRUE)

s3t2off <- stack(b3off,b2off) 
clusterR(s3t2off, overlay, args=list(fun=nd_fn),progress='text',filename="nd3t2off.tif", options=c("COMPRESS=DEFLATE"),datatype='INT2S', overwrite=TRUE)

s3t7off <- stack(b3off,b7off) 
clusterR(s3t7off, overlay, args=list(fun=nd_fn),progress='text',filename="nd3t7off.tif", options=c("COMPRESS=DEFLATE"),datatype='INT2S', overwrite=TRUE)

s4t5off <- stack(b4off,b5off) 
clusterR(s4t5off, overlay, args=list(fun=nd_fn),progress='text',filename="nd4t5off.tif", options=c("COMPRESS=DEFLATE"),datatype='INT2S', overwrite=TRUE)

s5t1off <- stack(b5off,b1off) 
clusterR(s5t1off, overlay, args=list(fun=nd_fn),progress='text',filename="nd5t1off.tif", options=c("COMPRESS=DEFLATE"),datatype='INT2S', overwrite=TRUE)

s5t4off <- stack(b5off,b4off)
clusterR(s5t4off, overlay, args=list(fun=nd_fn), progress='text', filename="nd5t4off.tif", options=c("COMPRESS=DEFLATE"),datatype='INT2S', overwrite=TRUE)

s7t3off <- stack(b7off,b3off) 
clusterR(s7t3off, overlay, args=list(fun=nd_fn),progress='text',filename="nd7t3off.tif", options=c("COMPRESS=DEFLATE"),datatype='INT2S', overwrite=TRUE)

s7t5off <- stack(b7off,b5off) 
clusterR(s7t5off, overlay, args=list(fun=nd_fn),progress='text',filename="nd7t5off.tif", options=c("COMPRESS=DEFLATE"),datatype='INT2S', overwrite=TRUE)

## Other normalized indices
# calcareous sediment index
calsed.off <- stack(b5off, b2off)
clusterR(calsed.off, overlay, args=list(fun=nd_fn),progress='text',filename="calsedoff.tif", options=c("COMPRESS=DEFLATE"),datatype='INT2S', overwrite=TRUE)

# ndvi - normalized difference vegitation index
ndvi.off <- stack(b4off, b3off)
clusterR(ndvi.off, overlay, args=list(fun=nd_fn),progress='text',filename="ndvioff.tif", options=c("COMPRESS=DEFLATE"),datatype='INT2S', overwrite=TRUE)


# msavi2 - modified soil adjusted vegitation index - for areas of sparce vegetation
msavi_fn <- function(bd1,bd2) {ind <- ((2*bd1+1)-(sqrt((2*bd1+1)^2-8*(bd1-bd2))))/2 
return(ind)
}
msavi.off <- stack(b4off, b3off)
clusterR(msavi.off, overlay, args=list(fun=msavi_fn),progress='text',filename="msavioff.tif", options=c("COMPRESS=DEFLATE"),datatype='INT2S', overwrite=TRUE)

# rock outcrop index
rock.off <- stack(b5off, b3off)
clusterR(rock.off, overlay, args=list(fun=nd_fn),progress='text',filename="rockoff.tif", options=c("COMPRESS=DEFLATE"),datatype='INT2S', overwrite=TRUE)



# end parallel cluster
endCluster()
gc()




