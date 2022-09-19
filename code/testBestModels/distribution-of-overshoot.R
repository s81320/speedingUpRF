# distribution of overshoot for models that were successful in test
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


################################################################
#### distribution of overshoot for models sucessful in test ####
################################################################
# which model is looked at : select by parameter X , line 75

doc.all <- data.frame(matrix(NA, nrow=0, ncol=5))
folder <- 'data/nursery03'
files <- list.files(folder, pattern='*.rda')#[1:2]

for(file in files){
  doc.f <- data.frame(matrix(0,nrow=0,ncol=ncol(doc.all)))
  
  # run loops over doc loaded from file
  load(paste(folder,file, sep='/')) # loads doc
  print(paste('loading',file,'at',Sys.time(), sep=' '))
  
  nReps <- length(doc)
  
  for(i in 1:nReps){
    #if(i%%10 ==0) print(paste(i/nReps , Sys.time()))
    
    DM <-  doc[[i]]$DM
    
    rg  <-  doc[[i]]$rg
    forest <-  rg$forest
    
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
    
    doc.f[i,] <- c(sz # size of (selected sub-) forest
    , chf$loRep
    , calcLogloss(subforest(forest, chf$forest), data.test) # logloss for (selecte sub-) forest
    , calcLogloss(subforest(forest, 1:sz), data.test) # regular small forest of the same size, reproducible
   # , calcLogloss(subforest(forest, sample(1:500,sz)), data.test) # regular small forest of the same size, not reproducible
    , calcLogloss(subforest(forest, order(LL)[1:sz]), data.test) # high performers of the same size
    )

  }
  doc.all <-  rbind(doc.all, doc.f)
  }
names(doc.all) <-  c('size','loRep','logloss.chf','logloss.bc.rs','logloss.bc.hp')

overshoot.bc.rs <- doc.all$logloss.chf - doc.all$logloss.bc.rs
par(mar=(c(4,4,1,1)+0.1))
hist(overshoot.bc.rs, breaks=100, main='',xlab='logloss overshoot', cex.axis=1.2, cex.lab=1.2)  
mean(overshoot.bc.rs)
shapiro.test(overshoot.bc.rs)
par(mar=(c(4,4,1,1)+0.1))
ggpubr::ggqqplot(overshoot.bc.rs) + 
  labs(x="", y="") + 
  theme(#axis.line=element_blank())
  axis.ticks=element_blank()
, axis.text.x=element_blank()
, axis.text.y=element_blank())
plot(density(overshoot.bc.rs), main='')

t.test(overshoot.bc.rs)

overshoot.bc.hp <- doc.all$logloss.chf - doc.all$logloss.bc.hp
hist(overshoot.bc.hp, breaks=100, main='',xlab='logloss overshoot', cex.axis=1.2, cex.lab=1.2)  
mean(overshoot.bc.hp)
plot(density(overshoot.bc.hp))

ggpubr::ggqqplot(overshoot.bc.hp) + 
  labs(x="", y="") + 
  theme(#axis.line=element_blank())
    axis.ticks=element_blank()
    , axis.text.x=element_blank()
    , axis.text.y=element_blank())

hist(doc.all$logloss.chf, breaks=100)
hist(doc.all$logloss.bc.rs, breaks=100)
hist(doc.all$logloss.bc.hp, breaks=100)
