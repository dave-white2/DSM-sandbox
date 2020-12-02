#12/2/20 correlation reduction, rfe



# Load and install packages
required.packages <- c("caret", "corrplot", "rgdal", "raster", "doParallel", "psych", "maptools", "dplyr", "sp", "snow", "snowfall", "plyr", "randomForest", "fmsb")
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(required.packages, require, character.only=T)
rm(required.packages, new.packages)


# Load Raster Data
# set working directory
setwd("C:/cov")

# read in raster layer names
rList=list.files(getwd(), pattern="tif$", full.names = FALSE)
rList

#create raster stack of covariates
rStack <- stack(rList)
names(rStack)

# create a dataframe
# for the following data reduction techniques, it is necessary to convert the raster stack into a data frame. This may cause problems with large raster stacks. my recommendation is to subset the raster stack by creating a regular sample of points. 
# If you need to create a regular sample due to a large set of covariates, that could be completed as follows:
df <- na.omit(as.data.frame(sampleRegular(rStack, na.rm=T, xy=F, size = 100))) #(length(rStack$calseddry)/2))))

plot(rStack$nd3t1dry)

length(rStack$calseddry)

gc()
#df <- na.omit(as.data.frame(rStack, na.rm=T, xy=F))


# filtering by correlation
# a correlation matrix is determined and covariates with high degree of correlation are returned

# create correlation matrix
corMat <- cor(df)

# visually examine the correlation matrix
corrplot(corMat, method = "circle")

# find high degree of correlation, cutoff is the threshold to set. If cutoff = 0.75 then covariates that are >= 75% correlated are removed
highCorr <- findCorrelation(corMat, cutoff = 0.85)

# remove covariates with high degree of correlation
df <- if(length(highCorr) > 0){
  df[, -highCorr]
} else {
  df
}

# subset rStack to crate a new raster stack of reduced covariates for use in other applications
r.stack <- subset(rStack, names(df))

## load in shapefile with points, if points are in a different folder change the dsn="." to the directory that it is in
shp.pts <-readOGR(dsn=".", layer="field_data_0206")


# the points and the raster stack need to be in the same projection
r.stack@crs # projection information for raster stack
shp.pts@proj4string # proj info for shapefile

all.equal(r.stack@crs, shp.pts@proj4string) #True - same projection

# if you need to convert the projection of the shape file, it would be done as follows:
# 
#  shp.pts <- spTransform(shp.pts, CRSobj =  r.stack@crs)


# r brings in a points shp file in as a SpatialPointsDataFrame (SPDF) object
# check the names of your SPDF
names(shp.pts)

# we only need the column that we are modeling remove all other columns
# remove columns 1 through 12 from shp.pts and rewrite it to shp.pts
shp.pts <- shp.pts[-c(1:12,14:18)]

# check the names of shp.pts to see if it worked you should only have the column with your modeling class
names(shp.pts)

## Plot to ensure alignment bw points and rasters
plot(r.stack$calsedwet)
plot(shp.pts, add=TRUE)


# this step gets the names of the rasters for use later
# convert raster stack to list of single raster layers
r.stack.list <- unstack(r.stack)
names(r.stack.list) <- names(r.stack)

# this following lines are the extract function this allows you to
# extract the covariate values at each data point collected

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

# create the training points by merging the extracted covariate values with the shape file
train.pts = merge(shp.pts, DF, by="ID")

#check the names to ensure merge
names(train.pts)

## Save points as a RDS file - this is a R file that contains an R object - save this for use later
#saveRDS(train.pts, file = "points.rds")

# clean junk out of memory
gc()

###################################################################
# this section gets the data ready for covariate reduction, where we select the covariates that are most usefule for modeling

# read in RDS file with points that include the extracted covariate values
#comp <- as.data.frame(readRDS("points.rds"))
comp <- as.data.frame(train.pts)

names(comp)

# remove the ID column and change the first column name to class
comp <- comp[,-c(1,43,42)] # removes the ID column
names(comp) # check the column names to make sure ID is removed
names(comp)[c(1)] <- c('Class') # change name of first col to class
names(comp) # check names of comp

# check the levels (classes) of comp class
levels(comp$Class)

# make sure Class is a factor

comp$Class <- as.factor(comp$Class)

# remove NA values from class
#comp <- comp[complete.cases(comp), ]

##################################################
# Recursive Feature Selection (RFE)
# this section is covariate reduction section


# This process is setup to run as a parallel
# in the make cluster function.

# Next, change the subsets
# to match the number of covariates that you have.

# check the number of covariates

length(comp) # number of covariates plus the class column

subsets <- c(1:(length(comp)-1))

# set seeds to get reporducable results when running the process in parallel
set.seed(238)
seeds <- vector(mode = "list", length=105)
for(i in 1:104) seeds[[i]] <- sample.int(1000, length(subsets) + 1)
seeds[[105]] <- sample.int(1000, 1)


# set up the rfe control
ctrl.RFE <- rfeControl(functions = rfFuncs,
                       method = "repeatedcv",
                       number = 15,
                       repeats = 5,
                       seeds = seeds, 
                       verbose = FALSE)

## highlight and run everything from c1 to stopCluster(c1) to run RFE

c1 <- makeCluster(detectCores()-1)
registerDoParallel(c1)
set.seed(9)
rf.RFE <- rfe(x = comp[,-1],
              y = comp$Class,
              sizes = subsets,
              rfeControl = ctrl.RFE,
              allowParallel = TRUE
)
stopCluster(c1)              

gc()

# 4.3  Look at the results
rf.RFE

#confusion matrix
rf.RFE$fit


#plotting rFE
plot(rf.RFE) # default plot is for Accuracy, but it can also be changed to Kappa
plot(rf.RFE, metric="Kappa", main='RFE Kappa')
plot(rf.RFE, metric="Accuracy", main='RFE Accuracy')

# see list of predictors
predictors(rf.RFE)

# the rfe function retuns the covariates with the highest accuracy for the rf model
# view the highest accuracy noted by the *

rf.RFE

# take the accuracy and subtract the accuracySD. look in the results of rf.RFE and find the accuracy that is > or = to this value. this is the number of covariates to use below

predictors(rf.RFE) # top number of covariates

# look at the top number of covariates that are equal to greater than the accuracy minus the accuracySD
predictors(rf.RFE)[1:5]

# assign this to a variable
a <- predictors(rf.RFE)[1:5]

# subset the covariate stack and the data frame by the selected covariates the variable a is your selected covariates

# subsed the raser stack to the selected covariates
r.stack.model <- subset(r.stack, a)

names(r.stack.model)

# subset the data frame points with the number of covariates selected
comp.sub <- (comp[,c("Class", a)])
names(comp.sub)


