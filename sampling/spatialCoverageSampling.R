# spatial coverage sampling

# load and install packages
required.packages <- c("raster", "rgdal", "spcosa")
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(required.packages, require, character.only=T)
rm(required.packages, new.packages)

setwd("C:/rwork/coweeta/cov")


# load in raster data
# read in raster layers names from modeling folder and create list
rlist=list.files(getwd(), pattern=".tif$", full.names = FALSE)
rlist

# create raster stack
rstk <- stack(rlist)
names(rstk)

# we only need one raster for this process
rstk <- rstk$chldem10m

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
writeOGR(samp.pts, dsn = "./clhs.pts", layer="SPcov135.pts", driver = "ESRI Shapefile", overwrite_layer = T)
