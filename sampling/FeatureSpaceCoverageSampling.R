# implementation of the feature space coverage sampling (FSCS) technique from:
# Ma, T., Brus, D. J., Zhu, A.-X., Zhang, L. and Scholten, T. (2020) Comparison of conditioned latin hypercube and feature space coverage sampling for predicting soil classes using simulation from soil maps. Geoderma, 370. https://doi.org/10.1016/j.geoderma.2020.114366


# load and install packages
required.packages <- c("raster", "rgdal", "LICORS", "parallel", "spcosa")
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

# create a subset raster stack of only selected covariates
# these would be selected from some sort of data reduction technique
names(rstk)
rstk <- subset(rstk, c("aspdeg", "chldem10m", "slpdeg", "tri"))
names(rstk)

# create a grid of points (n=gsize) to represent the raster stack
# - this can be done many ways sampleRegular is one method
gsize <- 10000 
samp.reg <- na.omit(sampleRegular(rstk, size = gsize, sp=T))

samp.reg.df <- na.omit(as.data.frame(samp.reg))
names(samp.reg.df)

# FSCS parameters
#nstart <- 100 # number of random starts
#sampNo <- 135 # number of clusters whose centroids represent sampling locations
#iter <- 1000 # number of iterations
# [,1:4] -columns representing covariates


# parallel implementation of FSCS
nw <- detectCores(logical = F)-1
cl <- makeCluster(nw)
clusterSetRNGStream(cl, iseed=1234)
set.seed(88)
samp.reg.df <- samp.reg.df
# Parallelize over the "nstart" argument
nstart <- 100
# Create vector of length "nw" where sum(nstartv) == nstart
nstartv <- rep(ceiling(nstart / nw), nw)
results <- clusterApply(cl, nstartv,
                        function(n,x) LICORS::kmeanspp(x, 135, nstart=n, iter.max=1000),
                        samp.reg.df[,1:4])
# Pick the best result
i <- sapply(results, function(result) result$tot.withinss)
result <- results[[which.min(i)]]
stopCluster(cl)

plot(result$inicial.centers)


#pull out the points
fscs <- as.data.frame(result$inicial.centers)
cov <- names(fscs)
fscs.pts <- merge(samp.reg.df, fscs, by = cov)

#projection information
prj <- rstk@crs
coords <- fscs.pts[c("x","y")]
#convert to sp
fscs.pts <- SpatialPointsDataFrame(coords = coords, data = fscs.pts, proj4string = prj)

plot(rstk$chldem10m)
points(fscs.pts, pch=19, col="red", cex=1.2)

#save pts to shpfile
writeOGR(fscs.pts, dsn = "./clhs.pts", layer="fscs135.pts", driver = "ESRI Shapefile", overwrite_layer = T)
