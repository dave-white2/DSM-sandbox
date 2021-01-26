# count of observations, (coobs) algorithm
# 
# # step 1 pixel to pixel comparison
# 1. calculate the multivariate distance to all other pixels in covariate stack
# 2. keep the max distance
# # step 2 pixel to point data comparison
# 1. calculate the multivariate distance from all pixels to pixel that represents the point (selected pixel)
# 2. standardize the distance
# 3. count the number of sites that are within a threshold for similarity to the selected pixel
#
# adapted from: https://bitbucket.org/brendo1001/clhc_sampling/src/master/additional/sample_tests_coobs.R
# Malone BP, Minansy B, Brungard C. 2019. Some methods to improve the utility of conditioned Latin hypercube sampling. PeerJ 7:e6451 https://doi.org/10.7717/peerj.6451

# load libraries
library(raster);library(sp); library(rgdal); library(clhs); library(snow); library(doParallel); library(parallel); library(clhs)

# set working directory
setwd("C:/rwork/clhs")

# import data
# point data
pts <- readOGR("ryanflat_pedons_all.shp")
# raster data
rlist = list.files(getwd(), pattern = "tif$", full.names = F)
# create raster stack
stack <- stack(rlist)
# aggregate for faster processing
stack <- aggregate(stack, fact = 6, fun=mean, na.rm=T)
# fact 6; 2.9 min
# fact 5; 6.3 min
# fact 3; 53.8 min
# no aggregation;long, long time

# coobs function - no editing necessary
coobs <- function(points, stack, prob){
  #data prep----
  att <- length(points)#checking for shpfile attributes
  pts <- points[, -(1:att)] #removing all attribute columns, only need location data
  ptDat <- as.data.frame(extract(stack, pts, na.rm=T))#extract covariate values to point data
  #set up raster data
  rDat <- as.data.frame(stack, xy=T) #convert to df
  rDat$cells <- cellFromXY(stack, cbind(rDat$x, rDat$y))
  names(rDat)
  rDat <- na.omit(rDat)
  # create new data frames
  rDatXY <- rDat[, c(1,2,(length(rDat)))] #extracting coordinate values
  rDat <- rDat[, c(3:(length(rDat)-1))] #keeping only covariate data
  names(rDatXY)
  names(rDat)
  #create function for progress bar
  f <- function(iterator){
    pb <- txtProgressBar(min = 1, max = iterator - 1, style = 3)
    count <- 0
    function(...) {
      count <<- count + length(list(...)) - 1
      setTxtProgressBar(pb, count)
      flush.console()
      c(...) # this can feed into .combine option of foreach
    }
  }
  iter <- nrow(rDat)
  covMat <- as.matrix(cov(rDat[,1:length(rDat)]))
  # make parallel cluster
  cpus = detectCores(logical = F)-1# the number of cores(not threads) in the cluster
  cl <- makeCluster(spec = cpus)
  registerDoParallel(cl)
  ### the actual function----
  oper1 <- foreach(i=1:iter, .combine = f(iter)) %dopar% {
    # pixel distance (raster)
    pix <- rDat[i, 1:length(rDat)]# pixel values
    pixDist <- mahalanobis(x = as.matrix(rDat[, 1:length(rDat)]), center = as.matrix(pix), cov = covMat) # calculate the distance of selected pixel to all other pixels
    minPix <- min(pixDist) # minimum distance (will alway be 0)
    maxPix <- quantile(pixDist, probs = prob) #maximum distance  ##the probs variable could change   ####### Hack to avoid outliers
    # point data distance
    datDist <- mahalanobis(x = as.matrix(ptDat[,1:length(ptDat)]), center = as.matrix(pix), cov = covMat)# calculate distance of observations to all other pixels
    datNdist<- (datDist-minPix)/(maxPix-minPix) # standardarise
    datNdist[datNdist > 1] <- 1 #if the datNdist is greater than 1 that means it datDist is greater than maxDist ##HACK
    datNdist <- 1- datNdist  # Higher values mean more similar
    #count how many obs are above a given threshold
    sum(datNdist >= prob)}
  #stop cluster 
  stopCluster(cl)
  rDatXY$sampleNOS <- oper1
  coob <- rasterFromXYZ(rDatXY[,c("x", "y", "sampleNOS")])
  return(coob)
}

# coobs function parameters
# coobs <- function(points, stack, prob)
# points -  your imported shapefile
#stack - your raster stack
#prob = 0.95 use as default#  #maximum distance; reducing would make more areas of similarity

#rCoobs <- coobs(points = pts, stack = stack, prob = 0.95)

# if you want to time how long it takes, run it the following way:
system.time({rCoobs <- coobs(points = pts, stack = stack, prob = 0.95)})

#plotting coobs raster
plot(rCoobs) 

#save the raster
#writeRaster(rCoobs, filename="sampleNos.tif", format="GTiff", overwrite=TRUE)

#using coobs as cost in clhs??
stack$cost <- rCoobs

#rs <- rasterToPoints(stack, spatial = T)
c50 <- clhs(stack, 100, cost='cost', iter=10000, progress=T, simple=F)

c50pts <- c50$sampled_data

plot(stack$cost)
points(c50pts, bg = 'red', pch=21)


# some notes
# 1. this is the original algorithm which employs the calculation of Mahalanobis distance, it could be possible to utilize other measures of similarity. 
# 2. aggregating the raster stacks beforehand is a way to cut down on calculation times, however this method should have more testing before implementation for the generation of points
# 3. alternatively, the authors suggest creating a sampling grid then of the covariate values and interpolate to the rest of the pixels, similar to the method I implemented and would still need more testing
# 4. adjusting the probability adjusts the 'similarity' threshold, higher probability is harder to find similar areas, lower probability includes more areas. to find the optimum range more testing is necessary, but a range of .8 to .99 should be sufficient. Potential rule of thumb, lower with less sample points??
