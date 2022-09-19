# depth of trees in simulation and in repetitions
# 16.9.2022

rm(list=ls())

source('code/source/distance-matrices.R')

# mean depth in trees in simulation (-> folder nursery02)
# and on full cleveland (-> folder nursery03)

# runt the following code twice:
# for sim=T and for sim=F
# in the end of this script: compare

sim <- T # check the simulations
if(sim){
  folder <-  'data/nursery02'
}else{
  folder <-  'data/nursery03' # repetitions on full Cleveland
}
files <- list.files(folder, pattern = '*.rda')#[1:3]
nFi <-  length(files)

doc1 <- list()

for(k in 1:nFi){ # for all files
  print(k)
  file <- files[k]
  load(paste(folder,file,sep='/')) # relevant data is in list doc

  N <- length(doc)
  doc2 <- matrix(0,N,14)  
  for(i in 1:N){
    rg <- doc[[i]]$rg
    DM <- doc[[i]]$DM
    # number of split nodes
  
    nTn <- tNodes(rg$forest,1) %>% length # number of terminal nodes
    fbt <- rt2fbt(rg$forest,1) # full binary tree
    nSn <- fbt %>% (function(x) which(x!=0)) %>% length # number of split nodes
    # assertthat::assert_that(nSn+1==nTn , msg='number of terminal nodes should be one more than number of split nodes.')
    mDissim <- Vectorize(function(j) mean(DM[[j]][upper.tri(DM[[j]])]))(1:4)
    min.max <- Vectorize(function(j) range(DM[[j]][upper.tri(DM[[j]])]))(1:4) %>% as.list %>% unlist
    doc2[i,] <-  c(nTn, log(x=length(fbt)+1, base=2) , mDissim , min.max)
  }
  doc1[[k]] <- data.frame(doc2)
  colnames(doc1[[k]]) <- c('number terminal nodes','height', paste('mean.',names(DM),sep=''), paste(c('min.','max.'),rep(names(DM),each=2),sep=''))
}

bind_rows(doc1) -> doc1

if(sim){
  doc.sim <- list(raw=doc1
                , meanValues= doc1 %>% apply(2,function(x) c(mean(x),sd(x)))
                , sim=sim)
}else{
  doc.fullCleve <- list(raw=doc1
                , meanValues= doc1 %>% apply(2,function(x) c(mean(x),sd(x)))
                , sim=sim)
}

doc.sim$meanValues ; doc.fullCleve$meanValues
doc.sim$meanValues / doc.fullCleve$meanValues


# result
# on simulated Cleveland and on original Cleveland
# trees have in the mean 22 terminal nodes (sd 2.5)
# and height 8
# pretty much the same, building the forest on simulations or the original Cleveland data

# a difference is in the standard deviation of the mean dissimilarity, it is around 10 for the d1, d2, Shannon Banks dissimilarities
# nothing to hint at why the d1 dissimilarity produces small Chipman based forests
# in tests, when built on the original Cleveland data


                