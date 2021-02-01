# project area stratification and generation of clhs points

# et up working environment
# load and install required packages ####
required.packages <- c("clhs","raster", "caret", "corrplot", "psych", "sp", "rgdal", "factoextra",
                       "doParallel", "foreach")
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(required.packages, require, character.only=T)
rm(required.packages, new.packages)

# set working directory for raster data
setwd("C:/rwork/coweeta/cov")

# load in raster data
# read in raster layers names from modeling folder and create list
rlist=list.files(getwd(), pattern=".tif$", full.names = FALSE)
rlist

# create raster stack
rstk <- stack(rlist)
names(rstk)


# create a subset raster stack of only selected covariates
names(rstk)
rstk.sub <- subset(rstk, c("aspdeg", "chldem10m", "slpdeg", "tri"))
names(rstk.sub)


# set working directory for shapefile
setwd("C:/rwork/coweeta")
# load in shapefile
poly <- readOGR("watershed.shp")


############################################################################################################
############################################################################################################


# stratify the project area utilizing only 4 selected covariates and generating 5 points per ----
# stratification area

# mask out raster stack by stratification area
# create a list to store masked raster stacks
out <- list()
# mask rstack by polygons in shapefile (creats a raster stack for each polygon)
# takes some time
for (i in 1:length(poly)){
  x <- crop(rstk.sub, poly[i,])
  out[[i]] <- mask(x, poly[i,])
}

rstk.spdf <- list()
for(n in 1:length(out)){
  rstk.spdf[[n]]<- rasterToPoints(out[[n]], spatial = T)
}

# set the number of clhs points you want generated
# this is the number of samples per stratification area
size = 5
# set the number of iterations for each clhs run
iter = 100000
# clhs points
# can take some time depending on the number of polygons
ss <- list()
for(o in 1: length(rstk.spdf)){
  ss[[o]] <- clhs(rstk.spdf[[o]], size = size, progress=T, iter=iter, simple=T)
}


clhs_pts <- bind(ss)# merges all spatialPointsDataFrames in list ss into one SPDF
plot(clhs_pts)
writeOGR(clhs_pts, dsn = "./clhs.pts", layer="pts135Strat", driver = "ESRI Shapefile", overwrite_layer = T)


############################################################################################################

###########################################################################################################

# what if we wanted to weight the points distribution by the area being sampled?
#

# inspecting the polygons for a field with area
poly
poly2 <- poly #making a copy of the polygons

# total number of samples we want to generate
totSaNo = 135


# manipulation to figure out how many samples in each stratification level
# this will produce more samples that the totSaNo due to rounding up to capture at least
# one point per area - this could use some further exploration

totAre <- sum(poly2$Area_ha) # total area
poly2$percent <- (poly2$Area_ha)/totAre # percent of area by each polygon
poly2$saNo <- (totSaNo) * (poly2$percent) # number of points per polygon based on area
poly2$saNo <- ceiling(poly2$saNo) # rounding so we have at least 1 point in every polygon


#clhs points
# the number of points is determined by poly2$saNo

# set the number of iterations for each clhs run
iter <- 100000
# run clhs on our stratified reduced raster stacks, utilizes the output, rstk.spdf
ss <- list()
for(o in 1: length(rstk.spdf)){
  s.size <- poly2$saNo[[o]]
  ss[[o]] <- clhs(rstk.spdf[[o]], size = s.size, progress=T, iter=iter, simple=T)
}


clhs_pts3 <- bind(ss) # merges all spatialPointsDataFrames in list ss into one SPDF
plot(clhs_pts3)
writeOGR(clhs_pts3, dsn = "./clhs.pts", layer="pts135StratWt", driver = "ESRI Shapefile", overwrite_layer = T)


#plot(clhs_pts, pch=18, col="red")
#points(clhs_pts2, pch=18, col="blue")
#points(clhs_pts3, pch=18, col="green")


###########################################################################################################
# generate 135 clhs points on entire area
rstk.sub.spdf <- rasterToPoints(rstk.sub, spatial = T)
ss<- clhs(rstk.sub.spdf, size = 135, progress=T, iter=iter, simple=T)
plot(ss)
writeOGR(ss, dsn = "./clhs.pts", layer="pts135", driver = "ESRI Shapefile", overwrite_layer = T)




