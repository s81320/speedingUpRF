# 10.8.2022
# code/chipmanSelection/basecases.R

# used in introduction
# tables tab:basecase:sim:1000 , tab:basecase:sim:100

# loops over simulation (files in folder nursery02) 
# to calculate loglosses for forests of different sizes


rm(list=ls())

library(dplyr)
library(ranger)
library(cluster)
library(diptest)
library(xtable)

load('data/data_SupMat4.rda') # loads the data sets Cleve, Hung, Swiss, VA
Swiss <- Swiss[,1:11]

source('code/source/prep.R') # calcLogloss (on forests) , calcLogloss2 (on predicted probabilities)
source('code/source/subforest.R') # subforest
source('code/source/chipman.R') # for calc_chipForest_2 , calc_LL_for_selection

basecases <- function(data, sizeSF, Nfiles=20, output=F){
  
  data.test <<- get(data) # data should be a variable name , character string 
  # get data from parent environment
  # is not passed to calc_LL_for_selection
  # is accessed from within the function (no good programming style!!)
  
  # to base the result on more bootstraps
  folder <- 'data/nursery02'
  files <- list.files(folder)[1:Nfiles] # restricting the number of simulations (in steps of 50)
  # dir(folder)
  
  LL.doc <- list()
  ct <- 1 # counter for list LL.doc
  for(file in files){
    if(output) print(paste(file, Sys.time()))
    # run loops over doc loaded from file
    load(paste(folder,file, sep='/')) # loads doc
    #print( parameter[p,])
    LL.doc.f <- data.frame(matrix(0, nrow=length(doc), ncol=length(sizeSF)))
    colnames(LL.doc.f) <- sizeSF
    # the following Vectorize is not the smart/fast way to do it
    # better calculate predictions for all trees (predict(forest, predict.all=T)) once
    # and then use predictions 1:x as predictions for the subforest
    # one might even build predictions for the forest of 1:(x+5) forest from predictions for the previous forest of 1:x trees and the new 5 trees
    for(i in 1:length(doc)){
      Vectorize(function(x) calcLogloss(forest=subforest(doc[[i]]$rg$forest,1:x), df=data.test))(sizeSF) ->
        LL.doc.f[i,]
    }
    LL.doc[[ct]] <-  LL.doc.f
    ct <- ct+1
    #if(output) print(paste('mean values in simulations in file', file,' : ', paste(round(LL.doc.f/length(doc),4),collapse=',')))
    if(output) print(apply(LL.doc.f,2,mean))
    #if(output) print(LL.doc)
  }
  
  LL.doc %>% bind_rows %>% data.frame -> doc.LL
  colnames(doc.LL) <-  sizeSF
  
  return(doc.LL)

}

# aggregate (calculate mean and sd) and
# generate LaTeX table
generate_table <- function(data){
  data %>% 
    #apply(2, function(x) c('mean'=median(x),'IQR'=IQR(x))) %>% 
    apply(2, function(x) c('mean'=mean(x),'sd'=sd(x))) %>% 
    t %>%
    xtable -> bc.xt
  digits(bc.xt) <- c(0,rep(4,2))
  #rownames(bc.xt) <- sizeSF # automatically uses row names from data (function argument)
  num.sim <- nrow(data)
  caption(bc.xt) <- paste('mean and standard deviation (over', num.sim,'simulations) for logloss for regular (small and default) forests of specified size' ,sep=' ')
  bc.xt %>% print
}

#### call function basecases for table of basecases:
#### loglosses for regular small forests (and the default forest)
#####################################################
sizeSF <- c(5,10,50,500)
bc100 <- basecases( data='Swiss',sizeSF, Nfiles=2 , output=F) # looping over 2 sim files, equals 100 simulations
bc1000 <- basecases( data='Swiss',sizeSF, Nfiles=20, output=F ) # looping over 20 sim files, equals 1000 simulations

generate_table(bc100)
generate_table(bc1000)

# small forests, sizes between 5 and 10
sizeSF <- c(5,6,7,8,9,10)
bc100 <- basecases( data='Swiss',sizeSF, Nfiles=2 , output=F) # looping over 2 sim files, equals 100 simulations
bc1000 <- basecases( data='Swiss',sizeSF, Nfiles=20, output=F ) # looping over 20 sim files, equals 1000 simulations

generate_table(bc100)
generate_table(bc1000)

# all sizes, nested, because we take the first trees (1:x) for a small forest , no sampling and random selection of trees
sizeSF <- c(seq(2,9,1),seq(10,500,5))
bc.all <-  basecases(data='Swiss',sizeSF, Nfiles=20 , output=F)
plot(sizeSF[-(1:4)],apply(bc.all,2,mean)[-(1:4)], type = 'l')

doc.LL <- list(sizeSF=sizeSF
               , logloss=bc.all
               , data.train='simulations of Cleveland'
               , data.test='Swiss'
               , N.sim=nrow(bc.all) # for the number of simulations 
               )
save(doc.LL, file='data/basecase/LL_nested_regular_small_08.2022_1000sim.rda')


