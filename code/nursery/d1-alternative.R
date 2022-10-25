# updated 25.10.2022 : new d1.alt.2 which is much faster than d1.alt.1

# new d1 formula ,tested for correctness
# with benchmarking
# subsubsection alternative definition sss:d1:altDef
# in chapter methods , section dissimilarities , subsection ss:d1

# microbenchmark at the end of the script
# parameters : number of trees (set directly in ranger call, num.trees) ,
#              nObs the number of observations needed to calculate d1

rm(list=ls())

source('code/source/distance-matrices.R')
source('code/source/plotTree.R')

load('data/CleveHungSwissVA.rda') # loads the data sets Cleve, Hung, Swiss, VA

library(ranger)
library(dplyr)
library(xtable)
library(microbenchmark)

set.seed(1)
rg <- ranger(CAD~.
             , data = Cleve[,1:11] 
             , num.trees = 500
             , replace = F 
             , mtry= 3 
             , importance = 'impurity'
             , probability = T 
             , min.node.size = 13 
)


#treeInfo(rg,4)
#plotTree1(rg,1)

nObs <- 100
data1 <- Cleve[1:nObs,]
pred.tn <- predict(rg, data=data1, type='terminalNodes')
ptnp <- pred.tn$predictions

createDMd1(forest=rg$forest, dft=data1)

################################################################################
# coded the definition for d1 quite verbatim 

# d1.def.1 working on the predictions (type terminalNodes)
d1.def.1 <- function(pred1, pred2){
  nObs <-  length(pred1)
  assertthat::assert_that(length(pred1)==length(pred2), msg='No can do. Need predictions of same length.')
  
  # masks / indicator matrix for observations in the same node for predictions 1
  I1 <- outer(1:nObs , 1:nObs , (function(i,j) pred1[i]==pred1[j])) 
  I1 <- I1[upper.tri(I1)]
  # ... and predictions 2
  I2 <- outer(1:nObs , 1:nObs , (function(i,j) pred2[i]==pred2[j]))
  I2 <- I2[upper.tri(I2)]
  
  # do predictions agree in putting observations in the same terminal node or not?
  # which pairs are treated the differently in predictions 1 and 2
  (I1!=I2) %>% factor(levels=c(TRUE,FALSE)) %>% table %>% prop.table %>% .['TRUE']
}

# d1.def calls d1.def.1
d1.def <-  function(rg , data){
  nObs <-  nrow(data)
  preds <- predict(rg, data, type='terminalNodes')$predictions
  outer(1:rg$num.trees 
        , 1:rg$num.trees 
        , Vectorize(function(i,j) d1.def.1(preds[,i],preds[,j])))
}

# returns the d1 dissimilarity matrix for all trees in rg
d1.def(rg, data1)

################################################################################
# my alternative definition

# pair_count counts the pairs given predictions of obs at terminal nodes
pair_count <- function(pred.tn){
  # predict the terminal nodes the obs in df fall to
  table(pred.tn) %>% unlist %>% unname %>% (function(x) x*(x-1)/2) %>% sum
}

pair_count(pred.tn$predictions[,2])

pair_count.idx <- function(pred.tn){
  #' returns indices for observations that form pairs in the terminal nodes
  nObs <-  length(pred.tn)
  I1 <- outer(1:nObs , 1:nObs , (function(i,j) pred.tn[i]==pred.tn[j])) 
  r1 <- which(I1==TRUE, arr.ind = TRUE)
  return(r1[(r1[,'row'] < r1[,'col']),])
}

pair_count.idx(pred.tn$predictions[,1])

f.alt <- function(pred1, pred2){
  nObs <-  length(pred1)
  assertthat::assert_that(length(pred1)==length(pred2), msg='No can do. Need predictions of same length.')
  I1 <- outer(1:nObs , 1:nObs , (function(i,j) pred1[i]==pred1[j])) 
  r1 <- which(I1==TRUE, arr.ind = TRUE)
  sameTn.idx.1 <- r1[(r1[,'row'] < r1[,'col']),]
  sameTn.idx.1 <-  matrix(sameTn.idx.1, ncol=2) # if sameTn.idx.1 has only 1 row it is interpreted as single numbers
  I2 <- outer(1:nObs , 1:nObs , (function(i,j) pred2[i]==pred2[j]))
  I2[sameTn.idx.1] %>% sum
}

f.alt(pred.tn$predictions[,4],pred.tn$predictions[,1])

# check this :
tri <- 5
f.alt(pred.tn$predictions[,tri],pred.tn$predictions[,tri])
pair_count(pred.tn$predictions[,tri])
# always equal

d1.alt.1 <-  function(ptnp, i,j){ # ptnp predict terminal nodes predictions
  nObs <-  length(ptnp[,i])
  ptnp.pc <- apply(ptnp, 2, pair_count)

  (pair_count(ptnp[,i]) 
    + pair_count(ptnp[,j]) 
    - 2*f.alt(ptnp[,i],ptnp[,j]))*2/(nObs*(nObs-1))
}

d1.alt.2 <-  function(ptnp, pc, i,j){ # ptnp predict terminal nodes predictions
  nObs <-  length(ptnp[,i]) # the same for all i ?!
  (pc[i]
    + pc[j] 
    - 2*f.alt(ptnp[,i],ptnp[,j]))*2/(nObs*(nObs-1))
}

# d1.alt.2 is much faster than d1.alt.1
pc <-  apply(ptnp,2,pair_count)
microbenchmark::microbenchmark(d1.alt.1(ptnp,1,2) , d1.alt.2(ptnp, pc, 1,2), times=2)

d1.alt <- function(rg , data){
  nObs <- nrow(data)
  ptnp <- predict(rg, data=data, type='terminalNodes')$predictions # P_red(... type = T_erminal N_odes)$ P_redictions
  pc <-  apply(ptnp,2,pair_count)
  
  A <- matrix(0,nrow=rg$num.trees, ncol=rg$num.trees)
  for(i in 2:rg$num.trees)
    for(j in 1:(i-1)){
      A[i,j] <- d1.alt.2(ptnp,pc,i,j) #  used to be d1.alt.1(ptnp,i,j) 
      A[j,i] <- A[i,j] # make it symmetric
    }
  return(A)
}

# this is far from optimal.
# f.alt builds I1 and calls pt.idx which also builds I1 , 
# => pair_count.idx should be integrated in f.alt
# I1 is built again and again when the same tree (index) is called, 
# in a forest of 100 trees it is built 99 times instead of once.
# and this happens for each of the 100 tree (in a forest of 100 trees)
################################################################################

# recall your parameters:
nObs
rg$num.trees

################################################################################
# compare runtimes

library(microbenchmark)

# check that they produce the same output
# check visually
d1.def(rg, data1)
d1.alt(rg, data1)
createDMd1(rg$forest, data1)

# check by code
assertthat::assert_that(all(d1.def(rg, data1)==createDMd1(rg$forest, data1)))
assertthat::assert_that(all(d1.def(rg, data1)==d1.alt(rg, data1)))


{ # the two implementations compare differently for 
  # num.trees 5, 50, 500 and number of observations 10, 100, 300
microbenchmark::microbenchmark(createDMd1(rg$forest, data1) 
              # , d1.def(rg, data1) 
               , d1.alt(rg, data1)
               , times=10
               , check='equal'
) -> mb.tab

mb.tab %>% print
ggplot2::autoplot(mb.tab) %>% print
mb.tab %>% summary %>% xtable(caption=paste(rg$num.trees,'trees,',nrow(data1),'observations'))
}


