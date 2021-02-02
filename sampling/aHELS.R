#  adapted hypercube evaluation of a legacy sample (aHELS)
#
# algorithm to allocate additional samples given existing samples
# basically it entails:
# 1. From the covariate data generate a matrix of quantiles (n x p) n = number of samples that are needed. p = number of available covariates
# 2. Using the quantile matrix create Hypercube matrices for both the existing point data and covariates.
# 3. Work out the densities of the elements in the hypercubes (count of obs/pixels within each quantile / total number of data or pixels)
# 4. Work out the ratio of sampling and grid densities... This will identify under and over sampling in the hypercube
# 5. To add additional samples:
#         1. rank the ratios from smallest to largest
#         2. workout the number of samples required to equalise quantile density of grids and sample data
#         3. Repeat step 5.2 until total number of additonal samples have been allocated.
#
#
# adapted from: https://bitbucket.org/brendo1001/clhc_sampling/src/master/additional/sample_tests_coobs.R
# Malone BP, Minansy B, Brungard C. 2019. Some methods to improve the utility of conditioned Latin hypercube sampling. PeerJ 7:e6451 https://doi.org/10.7717/peerj.6451


# libraries
library(raster);library(sp); library(rgdal)


# working directory
setwd("C:/rwork/clhs")

# input data
# read in raster layers and create list
rlist=list.files(getwd(), pattern="tif$", full.names = FALSE)
rlist

# create raster stack
r_stack <- stack(rlist)
names(r_stack)

# point data
pts <- readOGR("ryanflat_pedons_all.shp")


#### aHELS function

aHELS <- function(r_stack, pts, nosP){
  # raster data manipulation
  tempD <- as(r_stack, "SpatialPixelsDataFrame")
  tempD <- as.data.frame(tempD, xy=T, na.rm=T, cellnumbers = T)
  tempD$cellNos <- cellFromXY(r_stack, cbind(tempD$x, tempD$y))
  tempD <- na.omit(tempD)
  cols<- ncol(tempD)-3
  tempD <- tempD[, c((ncol(tempD)-2), (ncol(tempD)-1), ncol(tempD), 1:(ncol(tempD)-3))]
  
  
  #quantile matrix (of the covariate data)
  q.mat<- matrix(NA, nrow=(nosP+1), ncol= cols)
  j=1
  for (i in 4:ncol(tempD)){ #not the index start here
    #get a quantile matrix together of the covariates
    ran1<- max(tempD[,i]) - min(tempD[,i])
    step1<- ran1/nosP 
    q.mat[,j]<- seq(min(tempD[,i]), to = max(tempD[,i]), by =step1)
    j<- j+1}
  
  #covariate data hypercube
  cov.mat<- matrix(0, nrow=nosP, ncol=cols)
  n_iter <- nrow(tempD)
  # Initializes the progress bar
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = n_iter, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 100,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  
  for (i in 1:nrow(tempD)){ # the number of pixels
    cntj<- 1 
    for (j in 4:ncol(tempD)){ #for each column
      dd<- tempD[i,j]  
      for (k in 1:nosP){  #for each quantile
        kl<- q.mat[k, cntj] 
        ku<- q.mat[k+1, cntj] 
        if (dd >= kl & dd <= ku){cov.mat[k, cntj]<- cov.mat[k, cntj] + 1} 
      }
      cntj<- cntj+1
    }
    setTxtProgressBar(pb,i)
  }
  
  cov.mat[which(cov.mat==0)]<- 0.000001 # small number so we dont have to deal with zeros
  cov.mat # hypercube matrix
  covDens<- cov.mat/nrow(tempD)
  covDens # density matrix
  
  # point data manipulation
  att <- length(pts)#checking for shpfile attributes
  pts <- pts[, -(1:att)] #removing all attribute columns, only need location data
  
  # extract covariate values
  DSM_data<- extract(r_stack,pts, sp= T, method = "simple") 
  dat<- as.data.frame(DSM_data)
  dat<- dat[complete.cases(dat),]
  
  names(dat)[(length(dat)-1):(length(dat))]<- c("x", "y")
  dat$cellNos<- 0
  dat <- dat[, c((ncol(dat)-2), (ncol(dat)-1), ncol(dat), 1:(ncol(dat)-3))]
  
  #sample data hypercube (essentially the same script as for the grid data but just doing it on the sample data)
  h.mat<- matrix(0, nrow=nosP, ncol=cols)
  
  for (i in 1:nrow(dat)){ # the number of observations
    cntj<- 1 
    for (j in 4:ncol(dat)){ #for each column
      dd<- dat[i,j]  
      for (k in 1:nosP){  #for each quantile
        kl<- q.mat[k, cntj] 
        ku<- q.mat[k+1, cntj] 
        if (dd >= kl & dd <= ku){h.mat[k, cntj]<- h.mat[k, cntj] + 1}
      }
      cntj<- cntj+1
    }
  }
  
  # density
  h.mat[which(h.mat==0)]<- 0.000001  # small number so we dont have to deal with zeros
  h.mat 
  datDens<-h.mat/nrow(dat) # data density
  
  # selecting new samples
  rat<- datDens/covDens # ratio of data density and covariate density
  or<- order(rat) # rank the index where the biggest discrepancy is
  
  # indexes of quantiles that are not adequately sampled ie where rat is less than 1
  l1<- which(rat < 1, arr.ind = T)
  l1<- cbind(l1,which(rat < 1) )
  
  # What is the level of the greatest discrepancy? (This is important for comparing to later on when we have the additional sample)
  indy<- which(l1[,3]==or[1])
  rp<- l1[indy, 1]
  rc<- l1[indy, 2]
  rat[rp, rc]
  h.mat[rp, rc]
  cov.mat[rp, rc] # nuber of pixel (biggest discrepancy)
  datDens[rp, rc] # data density
  covDens[rp, rc] # covariate density
  
  # start from the highest discrepancy then work our way down
  upSamp<- nosP
  rpos<- 1
  base<- 3 # constant so that covariate columns are easily to select (realted to the column positions of the covariates in the tempD data frame)
  
  while (upSamp != 0){  # while the number of samples to allocate is greater than 0
    indy<- which(l1[,3]==or[rpos])
    rp<- l1[indy, 1]
    rc<- l1[indy, 2]
    
    ex<- floor(nrow(dat) * (datDens[rp,rc])) #existing count of samples within the selcted quantile
    eq<- ceiling(nrow(dat) * (covDens[rp,rc])) # number of samples needed to get to equal density between data and covariates
    sn<- eq-ex #number of samples needed
    if (upSamp < sn) {sn <- upSamp} # just so we dont over allocate
    
    
    #covariate selection
    covL<- q.mat[rp, rc]
    covU<- q.mat[rp+1, rc]
    subDat<- tempD[tempD[,(base+rc)] >= covL & tempD[,(base+rc)] <= covU,] # subset the covariates that meet the standard
    
    training <- sample( nrow(subDat), sn) #random number
    subDat2<- subDat[training,]
    
    #remove selcted samples from tempD so that repeated sampling does not occur (Is this necessary??)
    tempD<- tempD[!(tempD$cellNos %in% subDat2$cellNos), ]
    
    # Append new data to sampling dataframe
    dat<- rbind(dat,subDat2)
    
    #adjust the while params
    rpos<- rpos + 1 
    upSamp<- upSamp - sn
    print(sn)}
  
  # Check the sampling density with the addtional samples added
  #sample data hypercube
  h.mat<- matrix(0, nrow=nosP, ncol=cols)
  
  for (i in 1:nrow(dat)){ # the number of observations
    cntj<- 1 
    for (j in 4:ncol(dat)){ #for each column
      dd<- dat[i,j]  
      for (k in 1:nosP){  #for each quantile
        kl<- q.mat[k, cntj] 
        ku<- q.mat[k+1, cntj] 
        if (dd >= kl & dd <= ku){h.mat[k, cntj]<- h.mat[k, cntj] + 1}
      }
      cntj<- cntj+1
    }
  }
  
  #density
  h.mat[which(h.mat==0)]<- 0.000001
  h.mat
  datDens<-h.mat/nrow(dat)
  
  
  #### check
  rat<- datDens/covDens # ratio of data density and covariate density
  or<- order(rat) # rank the index where the biggest discrepancy is
  
  ## indexes of quantiles that are not adequately sampled
  l1<- which(rat < 1, arr.ind = T)
  l1<- cbind(l1,which(rat < 1) )
  
  # What the the level of the greatest discrepancy?
  indy<- which(l1[,3]==or[1])
  rp<- l1[indy, 1]
  rc<- l1[indy, 2]
  rat[rp, rc]
  h.mat[rp, rc]
  cov.mat[rp, rc]
  datDens[rp, rc]
  covDens[rp, rc]
  
  ## The following code does not have too much to do with the algorithm
  # Specify the different surveys (original and addtional)
  dat$survey<- NA
  str(dat)
  dat[dat$cellNos==0,"survey"]<- 1
  dat[dat$cellNos!=0,"survey"]<- 2
  
  ahels.pts <- subset(dat, survey == 2)
  
  ## Spatial points
  coordinates(ahels.pts) <- ~x + y
  str(ahels.pts)
  
  ## Coordinate reference systems
  proj4string(ahels.pts) <- pts@proj4string
  
  return(ahels.pts)
}

# aHELS function parameters
# aHELS <- function(r_stack, pts, nosP)
# r_stack - a raster stack of covariates
# pts - a SpatialPointsDataFrame containing field observations
# nosP - the number of sample points you want to genterate
ahels.pts <- aHELS(r_stack = r_stack, pts = pts, nosP = 100)


## Write point data to shapefile ; this shapefile will consist of both sampled points and new points
writeOGR(pts, ".", "hels100pts", "ESRI Shapefile")
