# script to test models with parameters found in hpo
# 18.5.2022

# around line 300 : generate effect plots for figure 8.1 labeled fig:test:lm:effects

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
# attribute data.set.name will be used in documentation 
# end of script, returned list of function f1

# function to do the tests
f3 <- function(testcases, folder='data/nursery03', num.files){
  #' test the given testcases on the forests and dissimilarity matrices in the given folder
  #' 
  #' @param testcases : dataframe giving methods, dissimilarities, purpose and more
  #' @param folder : where to get the forests and dissimilarity matrices from
  
  if(!is.null(num.files)){
    if(is.numeric(num.files)){
      files <- list.files(folder, pattern='*.rda')[1:num.files]
    }}else{
      files <- list.files(folder, pattern='*.rda')
    }
  
  doc.all <- data.frame(matrix(NA, nrow=0, ncol=9))
  
  for(file in files){
    doc.f <- data.frame(matrix(0,nrow=0,ncol=ncol(doc.all)))
    
    # run loops over doc loaded from file
    load(paste(folder,file, sep='/')) # loads doc
    print(paste('loading',file,'at',Sys.time(), sep=' '))
    
    nReps <- length(doc)
    
    for(i in 1:nReps){
      if(i%%10 ==0) print(paste(i/nReps , Sys.time()))
      
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
      
      # document test
      doctc <-  data.frame(matrix(NA,nrow=nrow(testcases),ncol=ncol(doc.all)))
      
      for(tc in 1:nrow(testcases)){
        
        method <- testcases[tc,'method']
        # use OOB observations for each tree to calculate the tree's logloss
        # prediction probabilities : pp
        
        dm <- DM[[testcases[tc,'dissim']]]
        
        parameter <- list('cutoff'=as.numeric(testcases[tc,'p.rep'])
                          , metric=testcases[tc,'dissim']
                          , 'sizeSF'= testcases[tc,'p.target.size']
                          , output=F
                          , selection='best')
        
        # a bit messy : some methods need the additional parameter selection 
        # some need oLL (ordered loglosses) but Meiner needs only LL 
        if(method=='chip1'){
          #        parameter$selection='best'
          grow_chipForest_1(dm = dm
                            , oLL= order(LL)
                            , parameter = parameter) -> chf 
        }
        
        if(method=='chip1_simplified'){
          #       parameter$selection='best'
          grow_chipForest_1_simplified(dm = dm
                                       , oLL= order(LL)
                                       , parameter = parameter) -> chf 
        }
        
        if(method =='meiner'){
          chf <- grow_meinForest(dm=dm
                                 , LL=LL
                                 , parameter=parameter)
        }
        
        if(method=='chip2'){
          chf <- grow_chipForest_2(dm=dm
                                   , oLL= order(LL)
                                   , parameter=parameter) # chipman 2 forest
        }
        
        # this method is totally suboptimal in this context, 
        # especially when the function is called for multiple sizes, like 50,10,9,8,7,6,5
        if(method=='hp'){
          chf <- list()
          chf$forest <- order(LL)[1:parameter$sizeSF]
          chf$loRep <- NA
        }
        
        sz <- length(chf$forest)
        
        # numerical data
        doctc[tc,c(1,4:8)] <- c( calcLogloss(forest, data.test) # full forest
                                 , testcases[tc,'p.rep'] # representation parameter
                                 , sz # size of (selected sub-) forest
                                 , chf$loRep
                                 , calcLogloss(subforest(forest, chf$forest), data.test) # logloss for (selecte sub-) forest
                                 , calcLogloss(subforest(forest, 1:sz), data.test) # logloss for regular small forest of same size as (selected sub-) forest
                                 #   , calcLogloss(subforest(forest,sample(1:forest$num.trees, sz)), data.test) # make the regular small forest random
        )
        
        # character , string
        doctc[tc,c(2,3,9)] <- testcases[tc,c('method','dissim','purpose')]
        
      }
      doc.f <-  rbind(doc.f, doctc)
    }
    doc.all <-  rbind(doc.all , doc.f)
  }
  names(doc.all) <- c('logloss.full.forest','method','dissimilarity','p.representation','size','alpha','logloss.chf','logloss.base.case','purpose')
  
  return(list('res'=doc.all, 'call'=list(testcases=testcases
                                         , folder=folder
                                         , num.files=num.files
                                         , data.train=attr(data.train,'data.set.name')
                                         , data.test=attr(data.test,'data.set.name')))
  )
}

#### run tests or load test previously run
# afterwards run evaluations and visualisations

#what_to_do <- 'run Chipman test'
#what_to_do <- 'run hp test'
what_to_do <- 'load tests'


################################################################
#### testing Chipman based forests with parameters from hpo ####
################################################################

if(what_to_do =='run Chipman test'){

  df.test <- data.frame(matrix(NA,nrow=0,ncol=5))
  names(df.test) <- c('method','purpose','dissim','p.rep','p.target.size')
  df.test
  # test Chipman 1 
  {
  df.prelim <- data.frame(matrix(NA,nrow=0,ncol=5))
  names(df.prelim) <- names(df.test)
  df.prelim[1,] <- c('chip1','mit50','d0',0,500)
  df.prelim[2,] <- c('chip1','freeLunch','d0',0.2,500)
  
  df.test <-  rbind(df.test, df.prelim)
}

  # test Chipman 2
  {
  df.prelim <- data.frame(matrix(NA,nrow=0,ncol=5))
  names(df.prelim) <- names(df.test)
  df.prelim[1,] <- c('chip2','freeLunch','d0',0,500)
  df.prelim[2,] <- c('chip2','mit50','d0',0,500)
  df.prelim[3,] <- c('chip2','mit5to10','d0',0.5,500)
  df.prelim[4,] <- c('chip2','mit5to10','d1',0.9,500)
  
  df.test <-  rbind(df.test, df.prelim)
}

  # test simplified Chipman 1, target size 50 and 5
  {
  df.prelim <- data.frame(matrix(NA,nrow=0,ncol=5))
  names(df.prelim) <- names(df.test)
  df.prelim[1,] <- c('chip1_simplified','mit50','d0',0,50)
  df.prelim[2,] <- c('chip1_simplified','mit50','d0',0.2,50)
  df.prelim[3,] <- c('chip1_simplified','mit50','d1',0.2,50)
  df.prelim[4,] <- c('chip1_simplified','mit50','sb',0,50)
  df.prelim[5,] <- c('chip1_simplified','mit50','sb',0.3,50)
  df.prelim[6,] <- c('chip1_simplified','mit5to10','d0',0,5)
  df.prelim[7,] <- c('chip1_simplified','mit5to10','d1',0.2,5)
  df.prelim[8,] <- c('chip1_simplified','mit5to10','d2',0.5,5)
  df.prelim[9,] <- c('chip1_simplified','mit5to10','sb',0.4,5)
  
  df.test <-  rbind(df.test, df.prelim)

}

  # test Meiner
  {
  df.prelim <- data.frame(matrix(NA,nrow=0,ncol=5))
  names(df.prelim) <- names(df.test)
  df.prelim[1,] <- c('meiner','freeLunch','d0',0,500)
  df.prelim[2,] <- c('meiner','freeLunch','d1',0.1,500)
  df.prelim[3,] <- c('meiner','freeLunch','sb',0.1,500)
  df.prelim[4,] <- c('meiner','mit50','d0',0,500)
  df.prelim[5,] <- c('meiner','mit50','d1',0.2,500)
  df.prelim[6,] <- c('meiner','mit50','sb',0.15,500)
  df.prelim[7,] <- c('meiner','mit5to10','d0',0.5,500)
  df.prelim[8,] <- c('meiner','mit5to10','d1',0.8,500)
  df.prelim[9,] <- c('meiner','mit5to10','d2',0.7,500)
  df.prelim[10,] <- c('meiner','mit5to10','sb',0.6,500)
  
  df.test <-  rbind(df.test, df.prelim)
}

  df.test$p.rep <- as.numeric(df.test$p.rep)
  df.test$p.target.size <- as.numeric(df.test$p.target.size) 
  #res <- f3(testcases=df.test[df.test$dissim=='d1',],num.files = 1)
  #res <- f3(testcases=df.test[df.test$method=='meiner',],num.files = 2)
  #res <- f3(testcases=df.test[df.test$purpose=='mit5to10',],num.files = 2)
  res <- f3(testcases=df.test,num.files = 20)
  
  #save(res, file='data/testNewMethods/test-all-methods-20files.rda')
}


###################################################
#### testing forests of high performing trees  ####
###################################################
# inside the if statement there is an option for which sizes of hp forests to test: 
# sizes 5,6,7,...49,50 (for X=1) or 5,7,9,..,499 (for other values of X)

if(what_to_do=='test hp'){
  # context: high performer and regular small forests as context for the Chipman based methods
 
  df.prelim <- data.frame(matrix(NA,nrow=0,ncol=5))
  names(df.prelim) <- names(df.test)
      
  X <- 1 # make your choice
  if(X==1){
    for(s in 5:50){
      df.prelim[s-4,] <-  c('hp',rep(NA,3),s)
      }
    }else{
      ct <- 1
      for(s in seq(5,500,2)){
        df.prelim[ct,] <-  c('hp',rep(NA,3),s)
        ct <- ct+1
      } # for
      } # else
  df.test.hp <- df.prelim
  df.prelim <- NULL
      
  df.test.hp$p.target.size <- as.numeric(df.test.hp$p.target.size)
  
  res <- f3(testcases=df.test.hp,num.files = 20)
  for(n in c('logloss.base.case','logloss.full.forest','logloss.chf','size')){
    res$res[,n] <-  as.numeric(res$res[,n])
    }
    
  #save(res, file='data/testNewMethods/test-hp-5to500-20files.rda')
  }
 
if(what_to_do=='load tests'){
  load('data/testNewMethods/test-all-methods-20files.rda') # loads res
}

#######################################################
#### now evaluate tests , generate plots / figures ####
#######################################################

r1 <- res$res

r1$overshoot.ff <- r1$logloss.chf - r1$logloss.full.forest
r1$overshoot.bc <- r1$logloss.chf - r1$logloss.base.case

#r1[r1$purpose=='mit5to10',] %>% 
  r1 %>%
    #group_by(purpose, method, dissimilarity, p.representation) %>% # the grouping decides about the order in the table we get in the end
    group_by(method, purpose, dissimilarity, p.representation) %>%
  summarize(m.size=mean(size), sd.size=sd(size), m.alpha=mean(alpha), m.chf=mean(logloss.chf), sd.chf=sd(logloss.chf), m.bc=mean(logloss.base.case), m.osh.bc=mean(overshoot.bc), m.osh.ff=mean(overshoot.ff)) -> tb1
View(tb1)

tb1 %>% 
  .[,names(tb1)[-7]] %>% # exclude alpha
  xtable(caption='test', label='tab:test:all',digits = c(rep(2,5),rep(1,2),rep(4,5)))

# values for purpose will be used in plot
# we change the values freeLunch -> avoid trade-off (size success trade-off) 
# mit50 -> mitigate50 (mitigate size success trade-off) , similar for mit5to10
{
  tb1[tb1$purpose=='freeLunch','purpose'] <-  'avoid trade-off'
  tb1[tb1$purpose=='mit50','purpose'] <-  'mitigate50'
  tb1[tb1$purpose=='mit5to10','purpose'] <-  'mitigate5to10'
}

# build linear models on aggregated data (tb1)
{
  # new variable overshoot, depending on purpose
  # for purpose 'avoid trade-off' it is the overshoot over the default forest's logloss,
  # for mitigation for size 50 and 5 to 10 it is the overshoot over the regular small forest of the same size
  tb1$m.osh <- ifelse(tb1$purpose=='avoid trade-off', tb1$m.osh.ff, tb1$m.osh.bc)  

  lm(formula='m.osh~method+purpose+dissimilarity', data=tb1) -> lm1 
  lm(formula='m.osh~method+purpose+dissimilarity+m.size', data=tb1) -> lm2 
  lm(formula='m.osh~method+purpose+dissimilarity+p.representation', data=tb1) -> lm3
  lm(formula='m.osh~method+purpose+dissimilarity+m.alpha', data=tb1) -> lm4
  lm1 %>% summary
  lm3$residuals %>% plot(type='b', main='residuals in linear model')

  # adding variables is not worth it
  print('anova with additional variable mean size')
  anova(lm1,lm2) %>% print
  print('anova with additional variable mean parameter of representation')
  anova(lm1,lm3) %>% print
  print('anova with additional variable mean level of representation')
  anova(lm1,lm4) %>% print

  # the following generates effect plots for figure 8.1 labeled fig:test:lm:effects
  # linear model on the aggregated data
  # focusing on data for a selected purpose (few observations!)
  effects::allEffects(lm1) %>% plot # 
  par(mar=c(4,4,1,1)+0.1)
  ylim=c(-0.05,0.02)
  effects::effect('method', lm1) %>% plot(ylab='overshoot', main='', ylim=ylim)
  effects::effect('dissimilarity', lm1) %>% plot(ylab='overshoot', main='', ylim=ylim)
  effects::effect('purpose', lm1) %>% plot(ylab='overshoot', main='', ylim=ylim)
}

# check if quadratic additional predictor is worth adding :
lm(formula='m.osh~method+purpose+dissimilarity+poly(m.size,2)', data=tb1) -> lm2 
lm(formula='m.osh~method+purpose+dissimilarity+poly(p.representation,2)', data=tb1) -> lm3
lm(formula='m.osh~method+purpose+dissimilarity+poly(m.alpha,2)', data=tb1) -> lm4
lm4 %>% summary
lm3$residuals %>% plot(type='b', main='residuals in linear model')
anova(lm1,lm4) %>% print

################################################################################
#### code below not in thesis ##################################################
################################################################################

{
  # tb2 <- tb1[tb1$purpose=='mitigate5to10',]
  #tb2 <- tb1[tb1$purpose=='avoid trade-off',]
  tb2 <- tb1[tb1$purpose=='mitigate50',]
  lm(formula='m.osh~method+dissimilarity', data=tb2) -> lm1 
  lm1 %>% summary %>% print
  lm1$residuals %>% plot(type='b', main='residuals')

  effects::allEffects(lm1) %>% plot # 
}

##### non-aggregated data ####
##############################


# new variable overshoot, depending on purpose
# for purpose freeLunch it is the overshoot over the default forest's logloss,
# for mitigation for size 50 and 5 to 10 it is the overshoot over the regular small forest of the same size
r1$overshoot <- ifelse(r1$purpose=='freeLunch', r1$overshoot.ff, r1$overshoot.bc)  
#r1$overshoot.bc <- NULL
#r1$overshoot.ff <-  NULL

# mixed models
lm(formula='overshoot~purpose+method+dissimilarity', data=r1) %>% summary
lm(formula='overshoot~purpose+method*dissimilarity', data=r1) %>% summary
lm(formula='overshoot~purpose*method+dissimilarity', data=r1) %>% summary
lm(formula='overshoot~purpose*dissimilarity+method', data=r1) %>% summary
lm(formula='overshoot~purpose*dissimilarity*method', data=r1) %>% summary
# they are all about the same
# of course there is a lot of variance that cannot be explained away (R squared small)
# the residual std error is reasonable
# and the p-value is very small

# check if the representation parameter gets a consistent level of representation // spoiler: it does, all is well !!
{lm(alpha~p.representation+dissimilarity, data=r1) %>% summary # Multiple R squared almost 1, p-value extremely small: great, res. std error : depends!
# alpha for Shannon-Banks has mean 29 and sd 4, so we can easily tolerate a res standard error of less than 1 : also great
lm(alpha~p.representation, data=r1[r1$dissimilarity=='sb',]) %>% summary
# the remaining dissimilarities have standard deviations 0.01 to 0.1, so a residual standard error of 0.04 is not really great
# we'd have to look at each dissimilarity separately
lm(alpha~p.representation+dissimilarity, data=r1[!r1$dissimilarity=='sb',]) %>% summary
}

# par(mar=c(4,4,2,1)+0.1)
# plot(lm1)

{
  d1 <- filter(r1,purpose=='mit50')
  lm(formula='overshoot~method+dissimilarity', data=d1) -> lm1
  lm1 %>% summary %>% print
  effects::allEffects(lm1) %>% plot
  #effects::Effect(mod = lm1,focal.predictors = 'dissimilarity') %>% plot
}
{ lm(formula='overshoot~method+purpose+dissimilarity', data=r1) -> lm1
  lm(formula='overshoot~method+purpose+dissimilarity+alpha', data=r1) -> lm2
  lm(formula='overshoot~method+purpose+dissimilarity+size', data=r1) -> lm3
  lm(formula='overshoot~method+purpose+dissimilarity*p.representation', data=r1) -> lm4
  anova(lm1,lm2)
  anova(lm1,lm3)
  anova(lm1,lm4)
  effects::allEffects(lm4) %>% plot(multiline=T)
}

summary(lm4)

row_selector <- which(r1$purpose=='mit5to10')
r1[row_selector ,c('size', 'overshoot.ff','overshoot.bc')] %>% apply(2, mean)
#r1[row_selector,c('purpose','overshoot.ff','overshoot.bc')] %>% View

# something with a regression line in r1
{
  row_selector <- which(r1$purpose=='mit5to10')
  lm1 <- lm(logloss.chf~ method+dissimilarity+logloss.full.forest+size , data=r1[row_selector,])
  lm1 %>% summary %>% print
#  lm2 <- lm(logloss.chf~ purpose+method+dissimilarity , data=r1)
#  lm2 %>% summary
  effects::allEffects(lm1) %>% plot
} # for mitigation to a size of 5 to 10, meiner is best

# drill into meiner, look at the dissimilarities
{
  row_selector <- which(r1$purpose=='mit5to10' & r1$method=='meiner')
  lm2 <- lm(logloss.chf~ dissimilarity , data=r1[row_selector,])
  lm2 %>% summary %>% print
  effects::allEffects(lm2) %>% plot
  lm3 <- lm(logloss.chf~ dissimilarity+size , data=r1[row_selector,])
  lm3 %>% summary %>% print
  effects::allEffects(lm3) %>% plot
  lm4 <-lm(size~dissimilarity, data=r1[row_selector,])
  lm4 %>% summary
  lm4 %>% allEffects %>% plot
} # dissimilarity Shannon-Banks is best (when we ignore size) with a model of p-value close to 1
# in linear model 3 when adding size, the p-value is great. And Shannon-Banks is worse, while d0 is best :-) 
# funny !!


# for mitigation of size-success-tradeoff for size around 50 trees
# linear model has p-value almost 1 , it's not worth anything!
{
  row_selector <- which(r1$purpose=='mitigate5to10')
  lm5 <- lm(overshoot.bc~ method+dissimilarity , data=r1[row_selector,])
  lm5 %>% summary %>% print
  #  lm2 <- lm(logloss.chf~ purpose+method+dissimilarity , data=r1)
  #  lm2 %>% summary
  effects::allEffects(lm5) %>% plot
}

# Meiner forest for a forest of size around 50
{
  row_selector <- which(r1$purpose=='mit5to10' & r1$method=='chip2')
  lm1 <- lm(overshoot.bc~ dissimilarity , data=r1[row_selector,])
  lm1 %>% summary %>% print
  #  lm2 <- lm(logloss.chf~ purpose+method+dissimilarity , data=r1)
  #  lm2 %>% summary
  effects::allEffects(lm1) %>% plot
}
# drill into meiner, look at the dissimilarities
{
  row_selector <- which(r1$purpose=='mit50' & r1$method=='meiner')
  lm1 <- lm(overshoot.bc~ dissimilarity , data=r1[row_selector,])
  lm1 %>% summary %>% print
  #  lm2 <- lm(logloss.chf~ purpose+method+dissimilarity , data=r1)
  #  lm2 %>% summary
  effects::allEffects(lm1) %>% plot
}


{r1 %>% 
    dplyr::filter(method=='meiner') %>%
    dplyr::filter(purpose=='mit5to10') %>% 
    lm(formula=logloss.chf~ logloss.full.forest + p.representation) -> lm3
    # lm(formula=overshoot.bc~ dissimilarity+ alpha ) -> lm3
    #lm(formula=size~ dissimilarity ) -> lm3
  
  lm3 %>% summary %>% print
  lm3 %>% 
    allEffects %>% 
    plot
}
