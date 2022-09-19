# cleveland02-v2/code/hyperparameters/pam_num_clusters.R
# Best number of clusters for dissimilarity matrices (based on each metric)
# needs to load dissimilarity matrices (does not create / calculate them)
# can load dissimilarity matrices based on Cleveland data or simulations thereof

# generate plots for figure 2.8 labeled fig:num.clus:1
# to check for the optimal number of clusters

rm(list=ls())

library(dplyr)
library(caret)
library(ranger)
library(effects)

withSim <- T
if(withSim){
  #load('data/nursery/nursery01_01_50x500.rda') # loads doc[[x]]$DM for x = 1:50
  load('data/nursery02/nursery02_01.rda')
  DM <-  doc[[1]]$DM
  }else{
    load('data/dms_for_a_forest_500_trees.rda') # loads DM
  }

metrices <-  c('d0','d1','d2','sb')

maxNumClus <-  50

# do clusterings and save silhouette widths
doc <- matrix(NA,maxNumClus-1,4)
for(m in 1:length(metrices)){
  metrices[m] %>% print
  dm <- DM[[m]]
  for(k in 2:maxNumClus){
    #print(k)
    pam.obj <-  cluster::pam(x=dm
                             , k=k
                             , diss=T
                             , nstart = 5
                             , keep.diss = F # faster
                             , keep.data = F #faster
                             )
    doc[k-1,m] <- pam.obj$silinfo$avg.width
  }
}

# plot silhouette width evolution for each dissimilarity
par(mar=c(2,2,1,1)+0.5)
for(m in 1:length(metrices)){
  plot(x=(1+(1:nrow(doc))) # clusters to 2 to 50
       , y=unlist(doc[,m])
       , type='l'
       #, main=paste('clustering based on metric' , metrices[[m]])
       #, xlab='number of clusters'
       , xlab=''
       #, ylab='average silhouette width'
       , ylab=''
       , ylim=c(0,1)
       , yaxt='n'
       , cex.axis=1.5
       )
  axis(side=2, at=seq(0,1,0.5),cex.axis=1.5)
  legend('topright'
         , legend=metrices[[m]]
         , cex=1.5
         , bty = 'n'
         )
  (which.max(doc[,m]) +1 ) %>% print
}


