#####
# Methods for generating sample points

# load and install packages
required.packages <- c("clhs", "raster", "rgdal", "spcosa", "doParallel")
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(required.packages, require, character.only=T)
rm(required.packages, new.packages)


#####################################################################################################################
setwd("C:/cov/sampling")


# load in raster data
# read in raster layers names from modeling folder and create list
rlist=list.files(getwd(), pattern=".tif$", full.names = FALSE)
rlist

# create raster stack
rstk <- stack(rlist)
names(rstk)

# we only need one raster for this process
rstk <- rstk$chldem10m




#####################################################################################################################
# regular grid sampling
# create a regular grid of points

regSample <- sampleRegular(rstk, size=150, sp=T)

#plotting
plot(rstk)
points(regSample, pch=19, col="red", cex=1.2)

# save points to shapefile
writeOGR(regSample, dsn = ".", layer="reg150.pts", driver = "ESRI Shapefile", overwrite_layer = T)




#####################################################################################################################
# random sample of grid
# create a random sample of points
rndmSample <- sampleRandom(rstk,  size=150, sp=T)

#plotting
plot(rstk)
points(rndmSample, pch=19, col="green", cex=1.2)

# save points to shapefile
writeOGR(rndmSample, dsn = ".", layer="rndm150.pts", driver = "ESRI Shapefile", overwrite_layer = T)




#####################################################################################################################
# spatial coverage sampling

# utilize the spatial coverage sampling from the spcosa package
# convert raster to data.frame
grd <- as.data.frame(rstk, xy=T)
grd <- na.omit(grd) #remove NAs
gridded(grd) <- ~ x * y # convert to grid format

#function
sc.grid <- stratify(grd, nStrata = 135, nTry=10)



#review
samp_pattern <- spsample(sc.grid)
plot(samp_pattern)

# pull out sample points and save to shp
samp.pts <- as(samp_pattern, "data.frame")
head(samp.pts)
coords<- samp.pts[c("x","y")]
prj <- rstk@crs
samp.pts <- SpatialPointsDataFrame(coords = coords, data = samp.pts, proj4string = prj)

#plotting
plot(rstk)
points(samp.pts, pch=19, col="blue", cex=1.2)

#save pts to shpfile
writeOGR(samp.pts, dsn = ".", layer="SPcov135.pts", driver = "ESRI Shapefile", overwrite_layer = T)





#####################################################################################################################
# conditioned latin hypercube sampling (cLHS)
#http://ncss-tech.github.io/soil-pit/sandbox/dave/clhs.html

#clhs utilizes covariate data to stratify the sampling points

# bring in raster data
setwd("C:/cov/sampling")


# load in raster data
# read in raster layers names from modeling folder and create list
rlist=list.files(getwd(), pattern=".tif$", full.names = FALSE)
rlist

# create raster stack
rstk <- stack(rlist)
names(rstk)

# convert the raster stack to points # Note for large datasets  use sample regular instead of converting the entire stack
stk.pts <- rasterToPoints(rstk, spatial=T)

c.pts <- clhs(stk.pts, size = 150, progress = T, iter=10000, simple = F, spatial=T)

# plot diagnostic plots
plot(c.pts, mode=c("obj", "box"))

plot(rstk$chldem10m)
points(c.pts$sampled_data, pch=19, col="blue", cex=1.2)

#save pts to shpfile
writeOGR(c.pts$sampled_data, dsn = ".", layer="clhs150.pts", driver = "ESRI Shapefile", overwrite_layer = T)


#####################################################################################################################
# conditioned latin hypercube sampling (cLHS) with cost
#http://ncss-tech.github.io/soil-pit/sandbox/dave/clhs.html

#clhs utilizes covariate data to stratify the sampling points

# bring in raster data
setwd("C:/cov/sampling")


# load in raster data
# read in raster layers names from modeling folder and create list
rlist=list.files(getwd(), pattern=".tif$", full.names = FALSE)
rlist

# create raster stack
rstk <- stack(rlist)
names(rstk)

# bring in cost raster
setwd("C:/cov/cost")
cost <- raster("cost.tif")

#add to raster stack
rstk$cost <- cost

# convert the raster stack to points # Note for large datasets  use sample regular instead of converting the entire stack
stk.pts <- rasterToPoints(rstk, spatial=T)

# note the implentation of cost in the clhs function
c.ptsCost <- clhs(stk.pts, size = 150, progress = T, iter=10000, cost='cost', simple = F, spatial=T)

# plot diagnostic plots
plot(c.ptsCost, mode=c("obj", "box"))

plot(rstk$chldem10m)
points(c.ptsCost$sampled_data, pch=19, col="blue", cex=1.2)

#save pts to shpfile
writeOGR(c.ptsCost$sampled_data, dsn = ".", layer="clhsCost.pts", driver = "ESRI Shapefile", overwrite_layer = T)





#####################################################################################################################
# spatial coverage sampling combined with clhs and cost

# utilize the spatial coverage sampling from the spcosa package
# convert raster to data.frame
grd <- as.data.frame(rstk, xy=T)
grd <- na.omit(grd) #remove NAs
gridded(grd) <- ~ x * y # convert to grid format

#function
sc.grid <- stratify(grd, nStrata = 300, nTry=10)

#review
samp_pattern <- spsample(sc.grid)
plot(samp_pattern)

# pull out sample points and save to shp
samp.pts <- as(samp_pattern, "data.frame")
head(samp.pts)
coords<- samp.pts[c("x","y")]
prj <- rstk@crs
samp.pts <- SpatialPointsDataFrame(coords = coords, data = samp.pts, proj4string = prj)

#plotting the spatial coverage points
plot(rstk$chldem10m)
points(samp.pts, pch=19, col="blue", cex=1.2)

#utilize spatial coverage points in the clhs function replacing the imput from the raster stack
# extract raster stack values to the spatial coverage points
samp.pts <- extract(rstk, samp.pts, na.rm = TRUE, df = TRUE, sp = TRUE)

# run clhs using the cost function
c.pts <- clhs(samp.pts, size=150, progress = T, iter = 10000, cost = 'cost', simple = F, spatial = T)

plot(rstk$chldem10m)
points(c.pts$sampled_data, pch=19, col="blue", cex=1.2)

#save pts to shpfile
writeOGR(c.pts$sampled_data, dsn = ".", layer="SPcovClhs.pts", driver = "ESRI Shapefile", overwrite_layer = T)


