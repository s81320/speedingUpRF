# code/hyperparameters/dip_unimodality.R
# used for table {tab:multimod:trees} table 2.8
# calculate the dip statistic

# how to plot results of a dip test is at the bottom, plot for convex and concave hulls of ecdf

rm(list=ls())

#install.packages('diptest')
library(diptest)
library(xtable)

# option A
#load('data/dms_for_a_forest_500_trees.rda') # loads DM

# option B
load('data/nursery02/nursery02_01.rda') # loads DM with 4 dissimilarity matrices, one for each dissimilarity

DM <- doc[[1]]$DM
names(DM)

for(metric in names(DM)){
  dm <-  DM[[metric]]
  calcDipPvalue <- function(i) dip.test(dm[i,-i])$p.value # uses dm from parent environment
  Vectorize(calcDipPvalue)(1:100) %>% hist(xlim=c(0,1), main=paste('dip-test on dissimlarity',metric,sep=' '), xlab='p-value')
}

# names(DM)
if(!exists('metrices')){
  if(is.null(names(DM))){
    metrices <-  c('d0','d1','d2','sb')
  }else{
    metrices <- names(DM)
  }
}

# generate data for table 2.8 , labeled as tab:multimod:trees
# set k.mds in loop to 2 , 50 , 100
par(mar=c(4,4,3,2)+0.1)
doc.dip <-  matrix(0,4,3)
for(i in 1:4){
  
  scaled.dm <- T
  if(scaled.dm){
    # p values for dip test on distance matrix for trees in simpler space / under mds
    k.mds <- 2
    DM[[i]] %>%
      cmdscale(eig = F, k = k.mds) %>% 
      dist %>% 
      as.matrix -> dm 
    stress <- sqrt( sum((DM[[i]]-dm)^2)/sum(dm^2) )
  }else{
    # p values for dip test on the dissimilarity matrices
    dm <- DM[[i]]
    stress <- NA
  }
  
  Vectorize(calcDipPvalue)(1:500) -> vdip # used dm defined above / in if statement
  vdip %>% 
    hist(xlim=c(0,1)
         , main=paste('histogram of p-values of dip statistic for',metrices[[i]])
    )
  doc.dip[i,] <- c(min(vdip), quantile(vdip,0.01), stress) 
  #ecdf(vdip)%>%plot(main=paste('empirical distribution function of p values of dip statistic\non all trees of forest with distances measured in',metrices[[i]]))
  #density(vdip)%>%plot(main=paste('density of p values of dip statistic on all trees of forest\ndissimilarity',metrices[[i]]))
}
colnames(doc.dip) <-  c('min p.value','1% quantile of p.values','stress')
doc.dip

"
#### visualize ####

mds<-function(D, xylim=FALSE ,main=NULL , col=NULL, pch=NULL){
  m <- cmdscale(D, eig = TRUE, k = 2)
  x <- m$points[, 1]
  y <- m$points[, 2]
  
  if(xylim){
    xlim<- c(-1,1)
    ylim<- c(-1,1)
  }else{
    xlim<-range(x)
    ylim<-range(y)
  }
  
  if(is.null(main)) main<-'mds'
  if(is.null(col)) col <- 1
  if(is.null(pch)) pch <- 1
  
  plot(x, y,  xlim = xlim, ylim=ylim , main=main, col=col, pch=pch)
  #text(x, y, pos = 4, labels =1:nrow(D))
  idx.min <- which.min(apply(D,2,sum))
  plot(x=x[c(1:3,idx.min)], 
       y=y[c(1:3,idx.min)],  
       xlim = xlim, 
       ylim=ylim , 
       main=main, 
       col=c(1,1,1,2), 
       pch=c(1,1,1,2))
}

for(m in 1:4){
  dm <- DM[[m]]
  
  doc <- data.frame(matrix(0,nrow(dm),2))
  
  for(i in 1:nrow(dm)){
    dt <-   dip.test(dm[i,-i])
    doc[i,] <-  c(dt$statistic, dt$p.value)
  }
  names(doc) <- c('dip','pval')
  
  # falling curve, not linear
  # plot(doc$dip , doc$p)
  
  plevel <-  0.01
  
  mds(dm 
    , main=paste('forest: mds for',metrices[[m]],'\nred : reject hypothesis of unimodality, p-level',plevel)
    , col=1+1*(doc$pval<plevel))
}

#### create a LaTex table ####

# create doc in long format to work well with dplyr from the tidyverse

doc <- data.frame(matrix(0,length(DM)*nrow(DM[[1]]),3))
ct <- 1 
for(m in 1:4){
  dm <- DM[[m]]
  
  for(i in 1:nrow(dm)){
    dt <-   dip.test(dm[i,-i])
    doc[(ct-1)*nrow(dm)+i,] <-  c(m,dt$statistic, dt$p.value)
  }
  ct <- ct+1
}
names(doc) <- c('metric','dip','pval')

doc$metric <- metrices[doc$metric]
doc %>% 
  group_by(metric) %>% 
  #filter(pval<0.05) %>%
  #filter(dip>0.1) %>%
  #summarise(mean(dip),mean(pval), length(which(pval<0.05)) , min(pval)) %>%
  summarise(length(which(pval<0.05)) , length(which(pval<0.01))) %>% 
  xtable -> xt

#digits(xt) <- 3
xt

dm <- DM[[2]]
dt <- dip(dm[i,-i], full.result=TRUE)
(dt$mj - dt$mn ) %>% hist(breaks=50)
dt$dip 
dip.test(dm)
plot(dt, main='')
"

