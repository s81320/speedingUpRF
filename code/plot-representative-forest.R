# interpreting the small forest by looking at its trees
# small forests of seizes between 5 and 10 can be interpreted by looking at its trees.

# we do this for the simplified Chipman 1 and Meiner forest under parameters that
# tested to be (in the mean) better than their basecases

# 7.9.2022

rm(list=ls())

library(ranger)
library(xtable)
library(dplyr)

load('data/data_SupMat4.rda') # loads the data sets Cleve, Hung, Swiss, VA

source('code/source/prep.R') # calcLogloss (on forests) , calcLogloss2 (on predicted probabilities)
source('code/source/subforest.R') # subforest
source('code/source/chipman.R') # for calc_chipForest_2 , calc_LL_for_selection
source('code/source/distance-matrices.R')
source('code/source/plotTree.R')


#### training data
{data.set.name <- 'Cleve'
  data.train <- get(data.set.name)[,1:11]
  attr(data.train,'data.set.name') <- data.set.name
  data.set.name <- NULL
}

#### test set
{data.set.name <- 'Hung'
  data.test <- get(data.set.name)[,1:11]
  attr(data.test,'data.set.name') <- data.set.name
  data.set.name <- NULL
}

load('data/nursery03/forests_and_DM_for_test_01.rda') # loads doc

i <- 1
DM <- doc[[i]]$DM
  
rg <- doc[[i]]$rg
forest <- rg$forest
  
predict(forest
        , data=data.train
        , predict.all = T)$predictions[,2,] -> pp # dim 303 x 500
  
# List of Loglosses : LL
lapply(1:forest$num.trees
         , function(t){ 
           OOB <- which(rg$inbag.counts[[t]]==0) 
           pp[OOB,t] %>% 
             calcLogloss2( df=data.train[OOB,] ) %>% 
             unlist
         }
       ) %>% 
  unlist -> LL
  
X <- 1
if(X=='1'){
  grow_chipForest_1_simplified(dm = DM[['d2']]
                                 , oLL= order(LL)
                                 , parameter = list('cutoff'=0.5
                                                    , 'sizeSF'=5
                                                    , 'selection'='best')) -> chf
  }else{
    grow_meinForest(dm=DM[['d0']]
                    , LL=LL 
                    , parameter=list('cutoff'=0.5, 'sizeSF'=500)) -> chf
    }

sz <-  length(chf$forest)
for(k in 1:sz){
  plotTree1(rg, chf$forest[k])
}

# variable importance in ranger
ranger::ranger(formula = 'CAD~.' 
               , data=Cleve[,1:11]
               , probability=T
               , importance = 'permutation') %>% 
  importance -> varImp

# which predictor variables are (most) important?
varImp[order(varImp)]
