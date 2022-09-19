# show differences in Chipman 1 forests with selection parameter best vs central
# chapter selecting subforests , hyperparameters,  hpo
# 9.8.2022

# chapter : new methods
# figure 7.2 in thesis

# compare and visualize d0 for best vs central
# for one / 1 forest , the one forest built on the first simulation of Cleveland

rm(list=ls())

library(ranger)
library(xtable)
library(caret)
library(dplyr)
library(diptest)
library(cluster)

load('data/data_SupMat4.rda') # loads the data sets Cleve, Hung, Swiss, VA

source('code/source/prep.R') # calcLogloss (on forests) , calcLogloss2 (on predicted probabilities)
source('code/source/subforest.R') # subforest
source('code/source/chipman.R') # for calc_chipForest_2 , calc_LL_for_selection
source('code/source/distance-matrices.R')
source('code/source/prep.R')


set.seed(321) # works well for comparing selection best vs central in Chipman 1
# using simulations we need the info on inbag observations per tree (why not use inbag of simulations?? because it is a different method)

load('data/sim1_with_inbag_info.rda') # loads doc2
load('data/nursery02/nursery02_01.rda') # loads doc
         
# training data is the 1st simulation of Cleveland
{data.train <- doc2$`bootstrapped training data`
attr(data.train,'data.set.name') <- 'Cleve_sim1'}
         
#### test set
{data.set.name <- 'Hung'
data.test <- get(data.set.name)[,1:11]
attr(data.test,'data.set.name') <- data.set.name
data.set.name <- NULL}
         
rg <- doc2$rg
forest <-  doc2$rg$forest
metric <- 'd0'
dm <- doc[[1]]$DM[[metric]] # get this from nursery simulation

# we use the dissimilarity matrix created for the forest in the nursery02 simulation with keep.inbag=F
# for the forest created on the same data with keep.inbag=T
# since we set the same seed in both cases, we get the same forests, and the same dissimilarity matrices.
# check by uncommenting the following assertion
# assertthat::assert_that((abs(createDMd0(rg$forest) - dm) %>% max)==0)

         
# logloss on Cleveland data not in simulation
OOB <-  Cleve[-doc2$resample,]

predict(forest
        , data=OOB
        , predict.all = T)$predictions[,2,] -> pp 
         
lapply(1:forest$num.trees
       , function(t) calcLogloss2(pp[,t], df=OOB) %>% unlist
       ) %>% 
  unlist -> LL

" # when working with the original Cleveland data set, no simulation (-> no uniform OOB set, need individual OOB / OOB for each tree individually)
# use OOB observations for each tree to calculate the tree's logloss
# prediction probabilities : pp
predict(forest
        , data=data.train
        , predict.all = T)$predictions[,2,] -> pp # dim 303 x 500

# List of Loglosses : LL
# logloss on individual trees OOB observations
lapply(1:forest$num.trees
       , function(t){ 
         OOB <- which(rg$inbag.counts[[t]]==0)
         pp[OOB,t] %>% 
           calcLogloss2( df=data.train[OOB,] ) %>% 
           unlist
       }
) %>% 
  unlist -> LL

"

parameter <- list('cutoff'=0.2) # cutoff for all forests, other parameters set individually

{ parameter$selection <- 'best'
  parameter$sizeSF <- 500
  cf1.best <- grow_chipForest_1(dm = dm
                                , oLL= order(LL)
                                , parameter = parameter
                                , output=T)
}

{ parameter$selection <- 'central'
  parameter$sizeSF <- 500
  cf1.central <- grow_chipForest_1(dm = dm
                    , oLL= order(LL)
                    , parameter = parameter
                    , output=T)
}

#################################################################
#### compare Chipman 1 forest with best and central tree per cluster

# top plot in figure 7.2
{par(mar=c(4,4,0,1)+0.2)
  # par(mar=c(4,4,2,1)+0.2) # when plotting a title
  plot(y=LL[order(LL)]
       , x=1:forest$num.trees
       , type='l'
       #, main='Logloss for all trees in forest, ordered'
       , xlim=c(0,500)
       , ylim=range(LL)
       , xlab='trees ordered by performance'
       , ylab='logloss (of individual / single tree)'
       , cex=1.5
       , cex.lab=1.5)
  
  oLL <-  order(LL)
  lapply(cf1.best$forest, function(i) which(i==oLL)) %>% unlist -> xXcf1.b
  points(x=xXcf1.b , y=LL[oLL[xXcf1.b]]+0.1, col='blue', pch='x' )
  #points(x=xXcf1.b , y=rep(0.5,length(xXcf1.b)), col='blue', pch='x' )
  
  lapply(cf1.central$forest, function(i) which(i==oLL)) %>% unlist -> xXcf1.c
  points(x=xXcf1.c , y=LL[oLL[xXcf1.c]]-0.1, col='blue')
  #points(x=xXcf1.c , y=rep(0.7, length(xXcf1.c)), col='blue')
  
#  lapply(cf1.central$denseR, function(i) which(i==oLL)) %>% unlist -> xXcf1.d
#  points(x=xXcf1.d , y=LL[oLL[xXcf1.d]], pch=0)
  
  legend('topleft', legend = c('Chipman 1, best tree per cluster','Chipman 1, central tree per cluster') 
         , col=c('blue','blue'), pch=c('x','o'), cex=1.5)
 # legend('bottomright' , legend=paste(metric, parameter$cutoff,sep=', '), cex=1.5)
}

# histograms in figure 7.2
# how each method / parameter setting for selection selects from the ordered trees
{hist(xXcf1.b, breaks=seq(0,500,50)) -> h.b.plt
hist(xXcf1.c, breaks=seq(0,500,50)) -> h.c.plt
ylim <- c(0, max(h.b.plt$counts , h.c.plt$counts))
par(mar=c(2,2,0,0)+0.4)
xXcf1.b %>% hist(breaks=seq(0,500,50)
                 , border='blue'
                 , xlim=c(0,500) 
                 , ylim=ylim
                 , xlab=''
                 , ylab=''
                # , xlab='bins of trees ordered by performance (bin size 50)'
                 , main = ''
                 , cex=1.5
                 #, cex.lab=1.2
                 , cex.axis=1.5)
legend('topright', legend='best', pch='x', col='blue', cex=1.5)
box()

xXcf1.c %>% hist(breaks= seq(0,500,50)
                       , border='blue'
                 , xlim=c(0,500) 
                 , ylim=ylim
                 , xlab=''
                 , ylab=''
                # , xlab='bins of trees ordered by performance (bin size 50)'
                 , main = ''
                 , cex=1.5
                 #, cex.lab=1.2
                 , cex.axis=1.5)
legend('topright', legend='central', pch='o', col='blue', cex=1.5)
box()
}

# data for bottom row table in figure 7.2
# mean logloss of single trees selected , logloss of these trees as a forest
# comparing Chipman 1 for selection parameter best vs central
LL[cf1.best$forest] %>% mean ; calcLogloss(forest=subforest(rg$forest, cf1.best$forest), df=data.test)
LL[cf1.central$forest] %>% mean ; calcLogloss(forest=subforest(rg$forest, cf1.central$forest), df=data.test)

length(cf1.best$forest) ; length(cf1.central$forest) ; cf1.best$size.dense.representing.sf

# clusters with only 1 tree => best is always central
cf1.best$hc.clus %>% table %>% (function(x) which(x==1)) %>% length
# clusters with 2 trees, for (around) half of them central should be best
cf1.best$hc.clus %>% table %>% (function(x) which(x==2)) %>% length

