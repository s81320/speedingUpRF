# This creates and saves forests and their dissimilarity matrices for testing
# always based on the full / original Cleveland data , no simulations
# needs keep.inbag=T to get individual tree OOB observations for sorting trees by performance
# 15.8.2022

rm(list=ls())

library(caret)
library(ranger)
library(dplyr)

metrices <- c('d0','d1','d2','sb')

load('data/data_SupMat4.rda')
# original Cleve is needed in createDMd1

source('code/source/distance-matrices.R')

doc <-  list()

# new seed for every new set of forests and dissimilarity matrices
# new seed for every file that is saved in the end
seed <- 20
set.seed(seed)

ct <- 1
nReps <-  50 # number of repetitions , number of forests saved at the end of this script
for(i in 1:nReps){
  #print(paste(i/nReps , Sys.time()))
  #if(i%%10 ==0) print(paste(i/nReps , Sys.time()))
  data.train <- Cleve[,1:11]
  
  rg <- ranger(CAD~.
               , data = data.train 
               , num.trees = 500
               , replace = F 
               , mtry= 3 # default
               , importance = 'impurity'
               , probability = T 
               , min.node.size = 13 
               , keep.inbag = T # needed for individual tree OOB
  )
  
  DM <- list()
  
  for(metric in metrices){
      # print(paste('    ' , metric , Sys.time()))
    DM[[metric]] <- createDM(forest=rg$forest , type=metric, dft=data.train)
  }
  
  doc[[i]] <-  list('data.train'='Cleve'
                    , 'rg'=rg 
                    , 'DM'=DM
                    , 'seed' = seed)
}

info <-  list('version'='1, training forests on original Cleveland data set and calculating the 4 dissimilarity matrices (d0, d1, d2, sb).'
              , 'seed'=seed
              , 'nReps'=nReps 
              , 'metrices'=metrices 
              , 'date'=Sys.time() 
              , 'created with script file'='code/nursery/02.R')
#save(info, doc, file='data/nursery03/forests_and_DM_for_test_***.rda')
