# show differences in (simplified) Chipman 1, Chipman 2 and Meiner forests
# for default forest built on 1st simulation of Cleveland
# review chapter in thesis

# generates figure 7.1 labeled fig:review:compare_models
# generates figure 7.2 labeled fig:review:compare_models_mds

# 9.8.2022

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

metric <- 'd1'

# generate or load a ranger random forest and specified representing subforests 
howTo <-  2
# 1 : generate from scratch, original Cleveland data, needs to calculate dissimilarity matrix (takes time, esp for d1)
# 2 : generate based on simulation of Cleveland, use dissimilarity matrix from nursery (faster)

set.seed(1234)
set.seed(321) # works well for comparing selection best vs central in Chipman 1
switch(EXPR=howTo, 
       {# generate from scratch, original Cleveland as training data for default forest
         # setting metric / dissimlarity before calculating the dissimilarity matrix
         #### training set
         {data.set.name <- 'Cleve'
         data.train <- get(data.set.name)[,1:11]
         attr(data.train,'data.set.name') <- data.set.name
         data.set.name <- NULL}
         # attribute data.set.name will be used in documentation (end of script)
         #### test set
         {data.set.name <- 'Hung'
         data.test <- get(data.set.name)[,1:11]
         attr(data.test,'data.set.name') <- data.set.name
         data.set.name <- NULL}
         # attribute data.set.name will be used in documentation 
         # end of script, returned list of function f1
         rg <- ranger(CAD~.
             , data = data.train
             , num.trees = 500
             , replace = F # T required for predictions ? Of course we want predictions!
             , mtry= 3 # default : 3
             , importance = 'impurity'
             , probability = T 
             , min.node.size = 13 
             , keep.inbag = T # needed to generate OOB obs as validation data
             )
         forest <- rg$forest
         dm <- createDM(forest=forest , type=metric , dft=data.train)
         
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
       } ,
       {# using simulations we need the info on inbag observations per tree (why not use inbag of simulations?? because it is a different method)
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
         dm <- doc[[1]]$DM[[metric]] # get this from nursery simulation
         
         # logloss on Cleveland data not in simulation
         OOB <-  Cleve[-doc2$resample,]
         predict(forest
                 , data=OOB
                 , predict.all = T)$predictions[,2,] -> pp 
         
         lapply(1:forest$num.trees
                , function(t){
                  calcLogloss2(pp[,t], df=OOB) %>% 
                    unlist}
         ) %>% 
           unlist -> LL
         } ,
      )

"# use OOB observations for each tree to calculate the tree's logloss
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
metric
parameter <- list('cutoff'=0.2)


{ # modify parameter (add selection, change sizeSF)
  parameter$selection <-  'best'
  parameter$sizeSF <-  50
  cf1.simp <- grow_chipForest_1_simplified(dm=dm
                                         , oLL=order(LL)
                                         , parameter=parameter
                                         , output=F)
  cf1.simp$forest
  cf1.simp$forest %>% length
}

{ parameter$sizeSF <- 500
  cf2 <- grow_chipForest_2(dm=dm
                         , oLL= order(LL)
                         , parameter=parameter) # chipman 2 forest
cf2$forest
cf2$forest %>% length 
}

{ parameter$sizeSF <- 500
  mf <- grow_meinForest(dm=dm
                      , LL=LL
                      , parameter=parameter
                      , output=F) # meiner forest
mf$forest
mf$forest %>% length
}

# sizes of representing subforests
length(cf1.simp$forest) ; length(cf2$forest) ; length(mf$forest)

# intersections, trees that appear in two subforests
base::intersect(cf1.simp$forest, cf2$forest) %>% length
base::intersect(mf$forest, cf2$forest) %>% length
base::intersect(cf1.simp$forest, mf$forest) %>% length

# trees that appear in all three subforests
base::intersect(cf1.simp$forest, cf2$forest) %>% 
  base::intersect(mf$forest) %>% 
  length

# plotting individual trees' loglosses for trees selected into representing sunforest
{par(mar=c(4,4,2,1)+0.1)
LL[cf1.simp$forest] %>% plot(type='p'
                             , col='green'
                             , main=''
                             , xlim=c(0,max(lengths(list(cf1.simp$forest, cf2$forest, mf$forest)))) 
                             , ylim=range(LL)
                             , xlab='tree'
                             , ylab='individual logloss')
LL[cf2$forest] %>% points(type='p', col='orange')
# LL[cf2$forest] %>% plot(type='p', main='Chipman 2' , xlim=c(0,500), ylim=range(LL))
oLLmf <- order(LL[mf$forest])
LL[mf$forest[oLLmf]] %>% points(type='p', col='red')
legend('topleft', legend=c('simplified Chipman 1','Chipman 2','Meiner'), col=c('green','orange','red'), pch='o')
}

#################################################################
#### compare Chipman 1 (best), simplified Chipman 1, Chipman 2 and Meiner forest 

# plot line for all trees in forest (500)
# plot dots for trees in representing forests 
{par(mar=c(4,4,0,1)+0.2)
  # par(mar=c(4,4,2,1)+0.2) # when plotting a title
plot(y=LL[order(LL)]
     , x=1:forest$num.trees
     , type='l'
     #, main='Logloss for all trees in forest, ordered'
     , xlim=c(0,500)
     , ylim=range(LL)
     , xlab='trees ordered by performance'
     , ylab='logloss (of individual / single tree)')

oLL <-  order(LL)

lapply(cf1.simp$forest, function(i) which(i==oLL)) %>% unlist -> xXcf1s
points(x=xXcf1s , y=LL[oLL[xXcf1s]]+0.1, col='green')

lapply(cf2$forest, function(i) which(i==oLL)) %>% unlist -> xXcf2
points(x=xXcf2 , y=LL[oLL[xXcf2]]-0.1, col='orange')

lapply(mf$forest, function(i) which(i==oLL)) %>% unlist -> xXmf
points(x=xXmf , y=LL[oLL[xXmf]]-0.2, col='red')

legend('topleft', legend = c('simplified Chipman 1','Chipman 2', 'Meiner') , col=c('green','orange','red'), pch=1)
legend('bottomright' , legend=paste(metric, parameter$cutoff,sep=', '))
}

# compare loglosses of default forest and representative subforests
calcLogloss(forest=rg$forest, df=data.test) # default
length(cf1.simp$forest) ; calcLogloss(forest=subforest(rg$forest, cf1.simp$forest), df=data.test)
length(cf2$forest) ; calcLogloss(forest=subforest(rg$forest, cf2$forest), df=data.test)
length(mf$forest) ; calcLogloss(forest=subforest(rg$forest, mf$forest), df=data.test)

#doc1 <- list()
#doc1[[5]] <- list(rg=rg , metric=metric , dm=dm , LL=LL , cf1=cf1 , cf2=cf2 , mf=mf)

#################################################################
#### compare Chipman 1 (best), simplified Chipman 1, Chipman 2 and Meiner forest 

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
       , cex.axis=1.5
       , cex.lab=1.5)
  
  oLL <-  order(LL)
  lapply(cf1.simp$forest, function(i) which(i==oLL)) %>% unlist -> xXcf1.s
  points(x=xXcf1.s , y=LL[oLL[xXcf1.s]]+0.1, col='green')
  #points(x=xXcf1.b , y=rep(0.5,length(xXcf1.b)), col='yellow', pch='x' )
  
  lapply(cf2$forest, function(i) which(i==oLL)) %>% unlist -> xXcf2
  points(x=xXcf2 , y=LL[oLL[xXcf2]]-0.1, col='orange')
  #points(x=xXcf1.c , y=rep(0.7, length(xXcf1.c)), col='yellow')
  
  lapply(mf$forest, function(i) which(i==oLL)) %>% unlist -> xXmf
  points(x=xXmf , y=LL[oLL[xXmf]]-0.2, col='red')
#  lapply(cf1.central$denseR, function(i) which(i==oLL)) %>% unlist -> xXcf1.d
#  points(x=xXcf1.d , y=LL[oLL[xXcf1.d]], pch=0)
  
  legend('topleft', legend = c('simplified Chipman 1','Chipman 2','Meiner') 
         , col=c('green','orange','red'), pch=1, cex=1.5)
 # legend('bottomright' , legend=paste(metric, parameter$cutoff,sep=', '), cex=1.5)
}

# barcharts
# how each method / parameter setting for selection selects from the ordered trees
{hist(xXcf1.s, breaks=seq(0,500,50)) -> hcf1s
hist(xXcf2, breaks=seq(0,500,50)) -> hcf2
hist(xXmf, breaks=seq(0,500,50)) -> hmf
ylim <- c(0, max(hcf1s$counts , hcf2$counts, hmf$counts))
par(mar=c(2,2,0,0)+0.7)

col <-  c('green','orange','red')
legend <-  c('simplified Chipman 1','Chipman 2','Meiner')
ct <- 1
for(xXdata in list(xXcf1.s, xXcf2, xXmf)){
  xXdata %>% hist(breaks=seq(0,500,50)
                  , border=col[ct]
                  , xlim=c(0,500) 
                  , ylim=ylim
                  , xlab=''
                  , ylab=''
                  # , xlab='bins of trees ordered by performance (bin size 50)'
                  , main = ''
                  , cex=2
                  #, cex.lab=1.2
                  , cex.axis=2)
  legend('topright', legend=legend[ct], pch='o', col=col[ct], cex=2)
  box()
  ct <- ct+1
}
}

# mean logloss of single trees selected , logloss of these trees as a forest
# comparing Chipman 1 for selection parameter best vs central
LL[cf1.simp$forest] %>% mean ; calcLogloss(forest=subforest(rg$forest, cf1.simp$forest), df=data.test)
LL[cf2$forest] %>% mean ; calcLogloss(forest=subforest(rg$forest, cf2$forest), df=data.test)
LL[mf$forest] %>% mean ; calcLogloss(forest=subforest(rg$forest, mf$forest), df=data.test)

length(cf1.simp$forest) ; length(cf2$forest) ; length(mf$forest) #cf1.best$size.dense.representing.sf


##########################################################################
#### compare models by plotting the subforests (mds)
#### probably not insightful as mds in k=2 dimensions has too much stress
#### and cannot represent original distances / dissimilarities

{mp <- cmdscale(d=dm, k=2) # for plotting we don't need eigenvalues, just the points
nT <- nrow(dm)
par(mar=c(1,1,1,1))
col <- c('green','orange','red')
legend <- c('simplified Chipman 1', 'Chipman 2', 'Meiner')
ct <- 1

plot(mp , col='black' , pch=1 , axes=F , xlab='' , ylab='', xlim= range(mp[,1]) , ylim=range(mp[,2]))
box(which = "plot", lty = "solid")
legend('bottomleft', legend=paste('tree in default forest', metric, sep=', ') , pch=1, cex=1.5)

for(f in list(cf1.simp$forest, cf2$forest, mf$forest)){
  plot(mp[f,] 
       , col=col[ct] 
       , pch=0 
       , axes=F 
       , xlab='' 
       , ylab='' 
       , xlim= range(mp[,1]) # needed to prevent automatic zooming in
       , ylim=range(mp[,2]) # needed to prevent automatic zooming in
  )
  legend('bottomleft', legend=legend[ct], pch=0, col=col[ct], cex=1.5)
  box(which = "plot", lty = "solid")
  ct <- ct+1
}
}

