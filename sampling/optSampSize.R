# optimum sample size when utilizing cLHS
# adapted from: https://bitbucket.org/brendo1001/clhc_sampling/src/master/additional/sample_tests_coobs.R
# Malone BP, Minansy B, Brungard C. 2019. Some methods to improve the utility of conditioned Latin hypercube sampling. PeerJ 7:e6451 https://doi.org/10.7717/peerj.6451


# set up
library(raster);library(rgdal);library(tripack); library(manipulate);library(clhs);library(entropy)

# working directory
setwd("C:/rwork/clhs")

# input data
# read in raster layers and create list
rlist=list.files(getwd(), pattern="tif$", full.names = FALSE)
rlist

# create raster stack
r_stack <- stack(rlist)
names(r_stack)


# function parameters
nb = 25 #Quantiles of the population # Number of bins # fewer bins less points needed
its<-2  # number internal iterations with each sample size number
cseq<-seq(5,500,50) # cLHC sample size
c.its <-100 #clHS iterations


# raster data setup
df <- as.data.frame(r_stack, na.rm=T)
ncols <- ncol(df)# number of covariates

#quantile matrix (of the covariate data)
q.mat<- matrix(NA, nrow=(nb+1), ncol= ncols)
j=1
for (i in 1:ncol(df)){ #note the index start here
  #get a quantile matrix together of the covariates
  ran1<- max(df[,i]) - min(df[,i])
  step1<- ran1/nb 
  q.mat[,j]<- seq(min(df[,i]), to = max(df[,i]), by =step1)
  j<- j+1}
#covariate data hypercube
cov.mat<- matrix(1, nrow=nb, ncol=ncols)
for (i in 1:nrow(df)){ # the number of pixels
  cntj<- 1 
  for (j in 1:ncol(df)){ #for each column
    dd<- df[i,j]  
    for (k in 1:nb){  #for each quantile
      kl<- q.mat[k, cntj] 
      ku<- q.mat[k+1, cntj] 
      if (dd >= kl & dd <= ku){cov.mat[k, cntj]<- cov.mat[k, cntj] + 1} 
    }
    cntj<- cntj+1
  }
}

#beginning of algorithm
mat.seq<- matrix(NA,ncol=2,nrow=length(cseq)) #empty matix for outputs
for (w in 1:length(cseq)){ # for every sample number configuration....
  s.size=cseq[w]  # sample size
  mat.f<- matrix(NA,ncol=2,nrow=its ) # placement for iteration outputs
  #internal loop
  for (j in 1:its){ #Note that this takes quite a while to run to completion
    repeat{
      ss <- clhs(df[,1:ncols], size = s.size, progress = T, iter = c.its)
      s.df<- df[ss,]
      if (sum(duplicated(s.df) | duplicated(s.df[nrow(s.df):1, ])[nrow(s.df):1]) < 2)
      {break}}
    ##Kullback-Leibler (KL) divergence
    h.mat<- matrix(1, nrow=nb, ncol=ncols)
    for (ii in 1:nrow(s.df)){ # the number of observations
      cntj<- 1 
      for (jj in 1:ncol(s.df)){ #for each column
        dd<- s.df[ii,jj]  
        for (kk in 1:nb){  #for each quantile
          kl<- q.mat[kk, cntj] 
          ku<- q.mat[kk+1, cntj] 
          if (dd >= kl & dd <= ku){h.mat[kk, cntj]<- h.mat[kk, cntj] + 1}
        }
        cntj<- cntj+1}}
    #h.mat 
    #Kullback-Leibler (KL) divergence
#    for(t in 1:ncols){
#      klo[t]<- KL.empirical(c(cov.mat[,t]), c(h.mat[,t]))
#      kloM <- mean(klo[1:t])
#      mat.f[j,1]<-kloM
#    }
    #Kullback-Leibler (KL) divergence
    klo.1<- KL.empirical(c(cov.mat[,1]), c(h.mat[,1])) #1
    klo.2<- KL.empirical(c(cov.mat[,2]), c(h.mat[,2])) #2
    klo.3<- KL.empirical(c(cov.mat[,3]), c(h.mat[,3])) #3
    klo.4<- KL.empirical(c(cov.mat[,4]), c(h.mat[,4])) #4
    klo.5<- KL.empirical(c(cov.mat[,5]), c(h.mat[,5])) #4
    klo.6<- KL.empirical(c(cov.mat[,6]), c(h.mat[,6])) #4
    klo<- mean(c(klo.1, klo.2,klo.3,klo.4,klo.5,klo.6))
    mat.f[j,1]<- klo  # value of 0 means no divergence
  } 
  #arrange outputs
  mat.seq[w,1]<-mean(mat.f[,1])
  mat.seq[w,2]<-sd(mat.f[,1])} ## END of LOOP

dat.seq<- as.data.frame(cbind(cseq,mat.seq))
names(dat.seq)<- c("samp_nos", "mean_KL","sd_KL")
##########################################################



#######################################################  
#plot some outputs
plot(cseq,mat.seq[,1],xlab="number of samples", ylab= "KL divergence")
plot(cseq,mat.seq[,2],xlab="number of samples", ylab= "standard deviation of percentage of total covariate variance of population account for in sample",main="Population and sample similarity")
#write.table(dat.seq, "Nav_datseq_clHC.txt", col.names=T, row.names=FALSE, sep=",")  # Save output to text file
# make plotting of the KL divergence
x<- dat.seq$samp_nos
y = 1 - (dat.seq$mean_KL)
plot(x, y, xlab="sample number", ylab= "CDF of 1-KL")  # Initial plot of the data



# Interactive plotting
start <- list()     # Initialize an empty list for the starting values
manipulate(
  {
    plot(x,y)
    a <- a; b <- b
    curve(a*x/(b+x), add = T)
    start <<- list(a=a, b=b)
  },
  a=slider(0.8, 1.2, step = .00001,  initial = .9),
  b=slider(5, 10, step = .00001, initial = 7)
)
fit1 <- nls(y~a*x/(b+x), start = start)
summary(fit1) 
lines(x, fitted(fit1), col="red")


#Apply fit
xx<- seq(1, 500,1)
lines(xx, predict(fit1,list(x=xx)))

jj<- predict(fit1,list(x=xx))
#normalized = 1- (jj-min(jj))/(max(jj)-min(jj))

x<- xx
y<- jj

plot(x, y, xlab="sample number", ylab= "cdf 1-KL", type="l", lwd=2)          # Initial plot of the data

x1<- c(-1, 500); y1<- c(0.95, 0.95)

lines(x1,y1, lwd=2, col="red")
x2 <- c(-1, 500); y1<- c(0.90, 0.90)
lines(x2,y1, lwd=2, col="red")
x.df <- data.frame(x)
y.df <- data.frame(y)
numbsamp <- cbind(x.df, y.df)


#############################################################################
