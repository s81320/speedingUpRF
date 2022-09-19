# introduction , methods , trees and random forests

# logloss performance for the single tree and 
# logloss for 3 individual trees vs the forest of these 3 trees 

# based on first 3 trees in first simulation

# table 3.2 in thesis , chapter methods , section trees and random forests

rm(list=ls())

file_name_script <-  'code/introduction/ensemble-magic.R'
list.files(path = 'code/introduction' , recursive = T)
assertthat::assert_that(file_name_script %in% list.files(recursive = T))

library(ranger)
library(xtable)
# library(caret)
library(dplyr)

load('data/data_SupMat4.rda') # loads the data sets Cleve, Hung, Swiss, VA

source('code/source/prep.R') # calcLogloss (on forests) , calcLogloss2 (on predicted probabilities)
source('code/source/subforest.R') # subforest
source('code/source/plotTree.R') # plotTree1 , plotTreeFb

load('data/sim1_with_inbag_info.rda') # loads doc2 and info

info
dft <- doc2$'bootstrapped training data' # data frame training
rg <- doc2$rg

data.unseen <- Swiss

rbind(
  # loglosses for single trees, each tree on its (individual) training data
  Inbag = c(calcLogloss(forest=subforest(rg$forest,1), dft[rg$inbag.counts[[1]]==1,]),
            calcLogloss(forest=subforest(rg$forest,2), dft[rg$inbag.counts[[2]]==1,]),
            calcLogloss(forest=subforest(rg$forest,3), dft[rg$inbag.counts[[3]]==1,]),
            NA , 
            NA
            ),
  # loglosses on training data (the forest's training data)
  training=c(calcLogloss(forest=subforest(rg$forest,1), dft),
             calcLogloss(forest=subforest(rg$forest,2), dft),
             calcLogloss(forest=subforest(rg$forest,3), dft),
             calcLogloss(forest=subforest(rg$forest,1:3), dft),
             calcLogloss(forest=rg$forest, dft)
             ) ,
  # loglosses on unseen data
  unseen=c(calcLogloss(forest=subforest(rg$forest,1), df = data.unseen),
           calcLogloss(forest=subforest(rg$forest,2), df = data.unseen),
           calcLogloss(forest=subforest(rg$forest,3), df = data.unseen),
           calcLogloss(forest=subforest(rg$forest,1:3), df = data.unseen),
           calcLogloss(forest=rg$forest, df = data.unseen)
           )
  ) %>% xtable 

# compare to best tree , CART

#install.packages('rpart')
library(rpart)
library(rpart.plot)

tri <- 1
data.train.inbag <- dft[which(rg$inbag.counts[[tri]]==1) ,] # Cleveland sim1 observations inbag for 1st tree
tree <- rpart::rpart(CAD ~ ., data = data.train.inbag) # same training data as 1st tree in 1st simulation

# predictions unpruned tree
predict(tree, data.train.inbag, type="prob")[,'Yes'] %>% calcLogloss2(df=data.train.inbag)
predict(tree, data.unseen, type="prob")[,'Yes'] %>% calcLogloss2(df=data.unseen)

# plot unpruned tree
rpart.plot::rpart.plot(tree)

# Step3: Prune the tree using the best cp.
bestcp <- tree$cptable[which.min(tree$cptable[,"xerror"]),"CP"]
tree.pruned <- rpart::prune(tree, cp = bestcp)

# plot pruned rpart tree
rpart.plot::rpart.plot(tree.pruned)

# predictions pruned tree
predict(tree.pruned, data.train.inbag, type="prob")[,'Yes'] %>% calcLogloss2(df=data.train.inbag)
predict(tree.pruned, data.unseen, type="prob")[,'Yes'] %>% calcLogloss2(df=data.unseen)





