# cLHS on points
# 
# cLHS can be implemented on a raster stack or a set of points, there just needs to be a dataframe passed to cLHS



# Load and install packages
required.packages <- c( "caret", "clhs", "e1071", "rgdal", "raster", "doParallel", "maptools", "dplyr", "sp", "snow", "snowfall")
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(required.packages, require, character.only=T)
rm(required.packages, new.packages)

# set directory that contains the covariates and point data
# *** pay attention to the / below
setwd("C:/cov")
getwd()


#Set up training data
# We need to bring the covariate data and site observations into R. Typically I keep all of the covariates and the point shapefile located in the same directory to keep things simple.

# Get the names of each covariate. All of mine are in the the .tif format. It is a good idea to work with consistent raster formats. My general prefernce for swithching between R, Saga, Arc, and others is the .tif format.

# read list of rasters
r.list=list.files(getwd(), pattern="tif$", full.names = FALSE) # reads file names from your working driectory that end in tif
r.list # view the list of rasters

#create raster stack of covariates
r.stack <- stack(r.list)

# check the names of the raster stack
names(r.stack)


## load in shapefile with training data points, if points are in a different folder change the dsn="." to the directory that it is in

#### IDF - only bring in the mountain slopes data. I would select all mountain slopes points and export into a new shapefile.
shp.pts <-readOGR(dsn=".", layer="field_data_0206")


# the points and the raster stack need to be in the same projection

# r brings in a points shp file in as a SpatialPointsDataFrame (SPDF) object
# check the names of your SPDF
names(shp.pts)

# we only need the column that we are modeling remove all other columns. In this case we are only interested in modeling the VegCLName. These were classes assigned during the field week, and roughly represent ecological sites.

# Remove the columns from the SPDF that we're not interested in at this time.
#shp.pts <- shp.pts[-c(1:12,14:18)]
 ### IDF - you will want to remove all columns in this case
  # try out the following code
shp.pts <- shp.pts[-c(1:length(names(shp.pts)))]


# check the names of shp.pts to see if it worked you should only have the column with your modeling class
names(shp.pts)
 ### IDF - you should have no columns at this point, bewlow we are going to extract the covariate values at each point

## Plot to ensure alignment bw points and rasters
# you can change the raster after the $ to reflect one of your rasters
plot(r.stack$rockdry)
plot(shp.pts, add=TRUE)


# this step gets the names of the rasters for use later
# convert raster stack to list of single raster layers
r.stack.list <- unstack(r.stack)
names(r.stack.list) <- names(r.stack)




# this following lines are the extract function this allows you to extract the covariate values at each data point collected

## Parallelized extract: (larger datasets)
sfInit(parallel=TRUE, cpus=parallel:::detectCores()-1)
sfLibrary(raster)
sfLibrary(rgdal)
# run parallelized 'extract' 
e.df <- sfSapply(r.stack.list, extract, y=shp.pts)
sfStop()
# clean memory
gc()
# now we need to assign the names of the covariates to the extracted values
DF <- as.data.frame(e.df)
names(DF) = tools::file_path_sans_ext(basename(names(r.stack.list)))
names(DF)
# head() looks at the top five rows of a data frame
head(DF)
# create ID field in both dataframe and shp
DF$ID <- seq.int(nrow(DF))
shp.pts$ID <- seq.int(nrow(shp.pts))

# create the points by merging the extracted covariate values with the shape file
pts = merge(shp.pts, DF, by="ID")

#check the names to ensure merge
names(pts)


## extract complete



#set up and run clhs

# apply clhs to points data frame
### IDF- size is the number of points that you want to keep
s.clhs <- clhs(pts, size = 250, progress = TRUE, iter = 10000, simple = FALSE)

# diagnostic plots
plot(s.clhs, mode = c("obj", "box"))

# extract indicies
subset.idx <- s.clhs$index_samples

# plot points to inspect and save shp file
# check visually:
par(mar = c(1,1,1,1))
plot(r.sagawi, axes=FALSE)
points(s[subset.idx, ], bg = 'red', pch=21)

# save cLHS points to shp
# change dsn to a working directory that you want to save to
# change layer name (file name)
writeOGR(s[subset.idx, ], dsn = 'C:/workspace2/clhs', layer = 'clhs_points_rgrd', driver = 'ESRI Shapefile', overwrite_layer = TRUE)

