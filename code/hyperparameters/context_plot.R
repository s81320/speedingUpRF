# visualize data generated in hyperparameters / ranger_num_trees-v3-nursery.R
# and from code / ChipmanSelection / 03_hpo.R

# used for figures 6.1 , 6.3 , 6.4
# fig:ch1:context:1 , fig:ch2:context:1 , fig:mei:context

library(dplyr)

# old code for old data
"{
  load('data/basecase/LL_nested_regular_forests_size_2_to_500.rda') # loads LL
  # what suits the data we have??
  apply(LL[1:100,],2,mean) -> LL.m # mean values only
  apply(LL,2,mean) -> LL.m # mean values only
}"

{
  # context for Chipman 1 : 100 simulations
  #load('data/basecase/LL_nested_regular_small_08.2022_100sim.rda') # loads doc.LL
  # context for Chipman 2 and Meiner: 1000 simulations
  load('data/basecase/LL_nested_regular_small_08.2022_1000sim.rda') # loads doc.LL
  LL <- data.frame(doc.LL$logloss)
  colnames(LL) <-  doc.LL$sizeSF # if need be
  LL %>% apply(2,mean) -> LL.m
}

folder <- 'data/chipman'
files <- list.files(folder) ; files
file <- files[[10]] ; file
load(paste(folder,file,sep='/')) # old hpo...rda files load et03 , 
# new ones load res1 with res1$et and additional info, for example on the method in res1$call$method
et <- res1$et # new hpo files
method <- res1$call$method # get the method from the data!

# not needed!
#for(cn in colnames(et)){
#  if(!cn=='selection')
#  et[,cn] <- as.numeric(et[,cn])
#}

# new methods in the context of regular random forests / the basecase

# rearrange data et to become et2
{
  colnames(et) # info on dissimilarity is in column names 
  # -> put dissim in a new column and stack columns for different dissimilarities
  # empty init
  matrix(nrow=0, ncol=5) %>%
    data.frame %>%
    setNames(c('logloss','I','size','dissim','cutoff')) -> et2

  metrices <- c('d0','d1','d2','sb')

  # fill et2:
  for(metric in metrices){
    et %>% 
      dplyr::select(ends_with(metric) | 'cutoff') %>%
      cbind('dissim'=metric) %>%
      setNames(c('logloss','I','size','cutoff','dissim')) %>%
      rbind(et2) -> et2
  }
}

# context plot for all sizes / full x-axis min.size:500
{ 
  min.size <- min(et2$size)
  colnames(LL) %>% as.integer %>% unlist %>% (function(x) x>min.size) %>% which -> idx
  ylim <- range(LL.m[idx], et2$logloss)
  #par(mar=c(4,4,2,1)+0.3)
  par(mar=c(4,4,1,1)+0.1)
  plot(colnames(LL[,idx]) %>% as.integer 
       , ylim=ylim
       , LL.m[idx] 
       , type='l'
       , xlab='size'
       , ylab='logloss'
    # , main=paste('success over size for', method,'forests\nin context of regular random forests', sep=' ')
       , cex.lab=1.2
       , cex.axis=1.2
    )
  points(et2$size , et2$logloss, col=as.factor(et2$dissim))
  legend('topright'
         , legend=c(paste(method, c('d0','d1','d2','sb'), sep=', '),'regular small forest')
       , pch=c(rep('o',4),'-')
       , col=c(1:4,1)
       , cex=1.2
  )
}

{
  par(mar=c(4,4,1,1)+0.1)
  size.lower.limit <- 5
  size.upper.limit <- 60

idx.zoom <-  which((as.integer(colnames(LL)) > size.lower.limit ) & (as.integer(colnames(LL)) < size.upper.limit ))
  LL.restricted <- LL[,idx.zoom] %>% apply(2,mean)
  ylim <- range(c(LL.restricted, et2[(et2$size>size.lower.limit)&(et2$size<size.upper.limit) ,'logloss']))

  LL.restricted  %>%
    plot(x=names(.) %>% as.integer
         , y=. 
         , ylim = ylim
         , type='l' 
         , xlab='size'
         , ylab='logloss'
         #, main=paste('success over size for', method,'forests\nin context of regular random forests', sep=' ')
         , cex.lab=1.2
         , cex.axis=1.2
       )

  idx <- which(et2$size > size.lower.limit) 
  text(et2[idx,'size']
     , et2[idx,'logloss' ]
     , labels = et2[idx,'cutoff' ] 
     #, cex=0.8
     , cex=1.2
     , col=as.integer(as.factor(et2[idx,'dissim'])))

  #legend('topright'
   #    , legend=c(paste(method, c('d0','d1','d2','sb'), sep=', '),'regular small forest', 'representation parameter')
    #   , pch=c(rep('o',4),'-','x')
     #  , col=c(1:4,1,1)
#       , cex=0.8
  #)
}

# context plot

method <- res1$call$method # get the method from the data!
size.lower.limit <- 4
size.upper.limit <- 6

ylim <- range(c(et2[points.selected,'logloss' ],LL.m[which((as.integer(colnames(LL))>size.lower.limit) & (as.integer(colnames(LL))<size.upper.limit ) )]))

LL[,which((as.integer(colnames(LL))>size.lower.limit) & (as.integer(colnames(LL))<size.upper.limit ) )] %>%
  apply(2,mean) %>%
  plot(x=names(.), y=. , 
       ylim=ylim,  type='l' , xlab='size', ylab='logloss')
points.selected <- (et2$size>size.lower.limit)&(et2$size<size.upper.limit)
points(et2[points.selected,'size'] 
       , et2[points.selected,'logloss' ]
       , col=as.factor(et2[points.selected,'dissim' ]))
legend('topright'
       , legend=c('regular RF', paste(method,metrices,sep=', '))
       , pch=c('-',rep('o',4))
       , col=c(1,1:4)
      # , cex=0.8
)


