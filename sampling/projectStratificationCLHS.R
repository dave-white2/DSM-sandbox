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


# stratify the project area utilizing only 4 selected covariates and generating 10 points per ----
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
size = 10
# set the number of iterations for each clhs run
iter = 10000
# clhs points
# can take some time depending on the number of polygons
ss <- list()
for(o in 1: length(rstk.spdf)){
  ss[[o]] <- clhs(rstk.spdf[[o]], size = size, progress=T, iter=iter, simple=T)
}


clhs_pts <- bind(ss)# merges all spatialPointsDataFrames in list ss into one SPDF
plot(clhs_pts)
writeOGR(clhs_pts, dsn = "./clhs.pts", layer="pts10cov4", driver = "ESRI Shapefile", overwrite_layer = T)


############################################################################################################
############################################################################################################


# what if you are interested in the covariate importance per stratification zone
# stratify by management zone, covariate reduction by zone, generate clhs points

# first remove covariates with nearZero variance
# convert raster stack to data frame
rstk.df <- as.data.frame(rstk, na.rm=T)

cl <- makeCluster(detectCores()-1) # makes a cluster using all but 1 cpu cores
registerDoParallel(cl)
# the actual zeroVar function, creates a vector containing which covariates should be removed
zeroVar <- nearZeroVar(rstk.df, foreach = TRUE, allowParallel = TRUE)
stopCluster(cl)

# removing unused items from memory
gc()

# check the zeroVar object
head(zeroVar)

# *note integer(0), means that there are no covariates to be removed

# remove covariates with zeroVar
rstk.df <- if(length(zeroVar) > 0){
  rstk.df[, -zeroVar]
} else {
  df
}

# remove the zerovar cov from raster stack
rstk.sub <- subset(rstk, names(rstk.df))


out <- list()
# mask rstack by polygons in shapefile (all covariates), creates a stack of all cov for
# each stratification zone
# takes some time
for (i in 1:length(poly)){
  x <- crop(rstk.sub, poly[i,])
  out[[i]] <- mask(x, poly[i,])
}

# convert list to data.frames
# time consuming step
out.df <- list()
for(j in 1:length(out)){
  out.df[[j]] <- as.data.frame(out[[j]], na.rm=T)
  }

#correlation reduction
out.hc <- list()
cutoff = 0.8
for(k in 1:length(out.df)){
  # filtering by correlation
  # Create correlation matrix
  out.hc[[k]] <- findCorrelation(cor(out.df[[k]]), cutoff = cutoff)
}


#subset to remove high corr vari
out.sub <- list()
for(l in 1:length(out.df)){
  out.sub[[l]] <- out.df[[l]][,-out.hc[[l]]]
}

# add ipca recuction when working properly
iPCA <- function(df, thresh = 95, ...){
  repeat{
    # create progress bar
    pb <- txtProgressBar(min = 0, max = 100, style = 3)
    # run pca on dataframe, can pass center and scale values
    pca.df <- prcomp(df, ...)
    sdev <- pca.df$sdev # extract the std dev
    evector <- pca.df$rotation # extract the eigen vectors
    evalue.df <- get_eigenvalue(pca.df) # extracts the eigenvalues and cumulative variance using the facto extra package
    evalue <- evalue.df$eigenvalue
    cum.var <- evalue.df$cumulative.variance.percent
    loadings <- as.data.frame((evector*evalue)/(sdev)) # calculate loadings per Jensen
    loadings$lf <- rowSums(abs(loadings))
    loadings <- loadings[order(-loadings$lf), ]
    loadings$cumVar <- cum.var
    idealthresh = thresh + 1
    score <- sum(tail(cum.var, 2))
    if(score >= (idealthresh + 100)){
      remove <- which(loadings$cumVar >= idealthresh)
      tokeep <- row.names(loadings[-c(remove),])
      df <- subset(df, select = tokeep)
      # progress bar for each iteration
      for(i in 1:100){
        Sys.sleep(0.1)
        # update progress bar
        setTxtProgressBar(pb, i)
      }
      close(pb)
    } 
    else if (score <= idealthresh + 100) break
  }
  return(df)
}

# The iPCA function 
out.ipca <- list()
for(m in 1:length(out.sub)){
  out.ipca[[m]] <- iPCA(out.sub[[m]], thresh = 95, center=T, scale=T)
  }

# subset raster stacks 
rstk.lst <- list()
for(o in 1:length(out.ipca)){
  rstk.lst[[o]] <- subset(out[[o]], names(out.ipca[[o]]))
}
# convert raster stacks to spatilaPointsDataFrames
rstk.spdf <- list()
for(p in 1:length(rstk.lst)){
  rstk.spdf[[p]]<- rasterToPoints(rstk.lst[[p]], spatial = T)
}

#clhs points
# this is the number of samples per stratification area
size = 10
# set the number of iterations for each clhs run
iter = 10000
ss <- list()
for(q in 1: length(rstk.spdf)){
  ss[[q]] <- clhs(rstk.spdf[[q]], size = size, progress=T, iter=iter, simple=T)
}



clhs_pts2 <- bind(ss)# merges all spatialPointsDataFrames in list ss into one SPDF
plot(clhs_pts2)

writeOGR(clhs_pts2, dsn = "./clhs.pts", layer="pts10covAllRed", driver = "ESRI Shapefile", overwrite_layer = T)




###########################################################################################################
###########################################################################################################



# what if we wanted to weight the points distribution by the area being sampled?
#

# inspecting the polygons for a field with area
poly
poly2 <- poly #making a copy of the polygons

# total number of samples we want to generate
totSaNo = 200


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
iter <- 10000
# run clhs on our stratified reduced raster stacks, utilizes the output, rstk.spdf
ss <- list()
for(o in 1: length(rstk.spdf)){
  s.size <- poly2$saNo[[o]]
  ss[[o]] <- clhs(rstk.spdf[[o]], size = s.size, progress=T, iter=iter, simple=T)
}


clhs_pts3 <- bind(ss) # merges all spatialPointsDataFrames in list ss into one SPDF
plot(clhs_pts3)
writeOGR(clhs_pts3, dsn = "./clhs.pts", layer="clhsSrtatWeight", driver = "ESRI Shapefile", overwrite_layer = T)



