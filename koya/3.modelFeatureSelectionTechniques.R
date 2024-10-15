# Modeling Feature Selection Techniques

# load and install package"s
required.packages <- c("Boruta","caret", "rgdal", "raster", "doParallel", "psych", "maptools", "dplyr", "sp", "snow", "snowfall", "plyr", "randomForest", "fmsb")
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(required.packages, require, character.only=T)
rm(required.packages, new.packages)




####################################################################################################################
# Data Prep
# set the working directory to where the covariate and pedon observation data are located
setwd('C:/cov')
getwd() # check to make sure the working directory is set

# read in the raster layers
rList=list.files(getwd(), pattern="tif$", full.names = FALSE)
rList

# create a raster stack
rStack <- stack(rList)
names(rStack)

# this step gets the names of the rasters for use later
# convert raster stack to list of single raster layers
rList
# remove the .tif
rList <- gsub(".tif", "", rList)
rList
names(rStack)
names(rStack) <- rList # assign names to the raster stack
names(rStack)

## load in shapefile with points, if points are in a different folder change the dsn="." to the directory that it is in
setwd('C:/cov')
shp.pts <- readOGR(dsn=".", layer="Sampled_Points")

# the points and the raster stack need to be in the same projection
crs(rStack)
shp.pts@proj4string # proj info for shapefile

all.equal(rStack@crs, shp.pts@proj4string) #True - same projection

# if you need to convert the projection of the shape file, it would be done as follows:
shp.pts <- spTransform(shp.pts, CRSobj =  rStack@crs)


# check the names of your SPDF
names(shp.pts)


# we only need the column that we are modeling remove all other columns
# remove columns 1 through 5 from shp.pts and rewrite it to shp.pts
shp.pts <- shp.pts[-c(1:5)]

# check the names of shp.pts to see if it worked you should only have the column with your modeling class
names(shp.pts)

# use the following line to chang your modeling col to Class for the script below
#names(shp.pts) <- "Class" ## changed name to Class for script below

# create plots to ensure alignment
plot(rStack$chldem10m)
plot(shp.pts, add=TRUE)


# this step gets the names of the rasters for use in the extract function below
# convert raster stack to list of single raster layers
r.stack.list <- unstack(rStack)
names(r.stack.list) <- names(rStack)
names(r.stack.list)


# this following lines are the extract function this allows you to
# extract the covariate values at each data point collected

## Parallelized extract: (larger datasets)
sfInit(parallel=TRUE, cpus=parallel:::detectCores()-4)
sfLibrary(raster)
sfLibrary(rgdal)
# run parallelized 'extract' 
e.df <- sfSapply(r.stack.list, extract, y=shp.pts)
sfStop()
# clean memory
gc()

# now we need to assign the names of the covariates to the extracted values
DF <- as.data.frame(e.df)
names(DF) = tools::file_path_sans_ext(basename(names(rStack)))
names(DF)
# head() looks at the top five rows of a data frame
head(DF)
# create ID field in both dataframe and shp
DF$ID <- seq.int(nrow(DF))
shp.pts$ID <- seq.int(nrow(shp.pts))

# create the training points by merging the extracted covariate values with the shape file
comp = merge(shp.pts, DF, by="ID")

#check the names to ensure merge
names(comp)

#Inspecting the spatial points data frame
str(comp)

# convert to a regular data frame
comp <- as.data.frame(comp, stringsAsFactors = T, spatial=T)
names(comp)
str(comp)

# remove any rows with NA values
comp <-comp[complete.cases(comp),]

# remove the ID col, the corrds cols and the spatial cols, these are from the shp file
# we only want the Class col and the covariate cols
names(comp)

# remove any unnecessary cols
comp <- comp[-c(1,48:50)] # these numbers should be changed depending on your dataset

names(comp) # you should be left with a col named Class and the rest are covariate data

comp$Class <- as.factor(comp$Class)# ensure Class is a factor
levels(comp$Class)

# clean junk out of memory
gc()




####################################################################################################################
# Recursive Feature Elimination (RFE)

# This process is setup to run as a parallel
# in the make cluster function.

# Set the number of subsets, this is related to the number of covariates
subsets <- c(1:length(comp)-1)

# set seeds to get reproducible results when running the process in parallel
set.seed(238)
seeds <- vector(mode = "list", length=112)
for(i in 1:111) seeds[[i]] <- sample.int(1000, length(subsets) + 1)
seeds[[112]] <- sample.int(1000, 1)

# set up the rfe control - adjustable parameters, see the caret package
ctrl.RFE <- rfeControl(functions = rfFuncs,
                       method = "repeatedcv",
                       number = 5,
                       repeats = 3,
                       seeds = seeds, 
                       verbose = FALSE)

## highlight and run everything from c1 to stopCluster(c1) to run RFE

c1 <- makeCluster(detectCores())
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

# Look at the results
rf.RFE

# look at confustion matrix with OOB
rf.RFE$fit

# this is the number of observations per class
summary(comp$Class)

# view the confusion matrix
rf.RFE$fit$confusion

# to save the convustion matrix run the following lines
#confusion.mat <- as.data.frame(rf.RFE$fit$confusion)
#write.table(confusion.mat, "RFEconfusion_matrix.txt", sep = "\t")

# look at the statistics for the results of RFE
rf.RFE$results

#plot(rf.RFE) # default plot is for Accuracy, but it can also be changed to Kappa
plot(rf.RFE, metric="Kappa", main='RFE Kappa')
plot(rf.RFE, metric="Accuracy", main='RFE Accuracy')

# see list of predictors, these are the selected covariates from RFE
predictors(rf.RFE)

# the rfe function retuns the covariates with the highest accuracy for the rf model
# view the highest accuracy noted by the *
rf.RFE

predictors(rf.RFE) # top number of covariates


# assign the top predictors to a variable to reduce our raster stack
a <- predictors(rf.RFE)

# subset the covariate stack and the data frame by the selected covariates the variable a is your selected covariates
# subset the raster stack to the selected covariates
rStackModel <- subset(rStack, a)

names(rStackModel)




#######################################################################################################################
# Boruta
# run the Boruta algorithm
fs_bor <- Boruta(y = comp$Class, x = comp[, -1], maxRuns = 35, doTrace = 1)

# plot variable importance and selected features
plot(fs_bor)

# plot evolution of the feature selection
plotImpHistory(fs_bor)

# extract the selected feature variables
vars <- getSelectedAttributes(fs_bor)

# view summary of the results
View(attStats(fs_bor))

# subset raster stack
rStackModel2 <- subset(rStack, vars)

names(rStackModel2)
