# 2.5.2022

# definitions of representation and diversity.
# figure 5.1 labeled fig:chip:representation:1
# figure 5.2 labeled fig:chip:diverse:1
# figure 5.4 labeled fig:chip:stress:1

# visualisation for representation and diversity

rm(list=ls())

library(MASS)
library(ranger)
library(dplyr)


plot_mds_forest <- function(D,m){
  x <- m$points[, 1]
  y <- m$points[, 2]
  
  idx.min <- which.min(apply(D,2,sum))
  
  main <- paste('mds of full forest', metric)
  col <-  rep(1, ncol(D))
  col[idx.min] <- 2
  pch <-  rep(1, ncol(D))
  pch[idx.min] <- 15
  plot(x=x, y=y, main=main, col=col, pch=pch)
  # text(x, y, pos = 4, labels =1:nrow(D))
}

# representation

plot_mds_representation <- function(D,m,qu.arg){
  x <- m$points[, 1]
  y <- m$points[, 2]
  
  idx.min <- which.min(apply(D,2,sum))
  
  # level of representation
  lor <- quantile(as.dist(D),qu.arg)
  # as.dist does not contain elements on the diagonal 
  # assertthat::assert_that((D %>% as.dist %>% length)==nrow(D)*(nrow(D)-1)/2)
  
  plot(density(D[idx.min,-idx.min])
     , main=paste('density for dissimilarity to most common tree', metric))
  abline(v=lor, col='grey') # line should end at density curve

  col <- rep(1,ncol(D))
  col[which(D[idx.min,]<lor)] <- 3
  col[idx.min] <- 2
  main <- paste('trees represented by central tree under the' 
              , metric
              , 'dissimilarity\nat level of representation'
              , round(lor,4) 
              , '(parameter ', qu.arg,')')
  pch <-  rep(1, ncol(D))
  pch[which(D[idx.min,]<lor)] <- 1
  pch[idx.min] <- 15

  par(mar=c(4,4,1,1)+0.1)
  plot(x=x
     , y=y
   #  , main=main
     , col=col
     , pch=pch
   , xlab='x mds'
   , ylab ='y mds'
   #, cex.lab=1.2
   )
  legend('bottomleft'
       , legend=c('central tree', 'tree represented by central tree', 'tree not represented by central tree')
       , col=c(2,3,1)
       , pch=c(15,1,1)
       #, cex = 1.2
       , bty='n'
       , box.lwd = 0
       )  

  # number of represented trees:in column TRUE
  (D[idx.min,]<lor) %>% table
  
  return(list('trees_represented'= which(D[idx.min,]<lor) 
              , 'central_tree'=idx.min
              , 'metric'=attr(D,'metric')))
}

stress <-  function(D,m){
(as.dist(D) - dist(m$points) )^2 %>% 
  sum %>% 
  (function(x) x/sum(as.dist(D)^2)) %>% 
  sqrt
}

####################################################
#### diversity
#### diverse representing subforest based on mds

plot_mds_div_rep_sf <-  function(D,m, qu.arg){
  
  x <- m$points[, 1]
  y <- m$points[, 2]
  
  idx.min <- which.min(apply(D,2,sum))
  
  # qu.arg is the parameter to set the level of diversity , applying it to the quantile function
  lod <- quantile(as.dist(D),qu.arg) # level of diversity

  drsf <- c(1) # diverse representing sub forest
  for(i in 2:ncol(D)){
    if(min(D[drsf,i])>lod){
      # i is not yet represented by drsf
      drsf <- c(drsf,i) # then add i to representing subforest
      # and drsf will still be diverse
    }
  }
  col <- rep('green',ncol(D)) # finally all trees wil be represented
  col[drsf] <- 'blue' # the representing trees
  pch <-  rep(1, ncol(D)) # crcle
  pch[drsf] <- 0 # 8 is the snowflake, 0 is the square
  main <- paste('diverse representing subforest', metric,'\n(level of diversity', round(lod,3),', parameter', qu.arg,')')
  plot(x=x, y=y
    # , main=main
     , col=col, pch=pch, xlab=' mds x', ylab='mds y')
  legend('bottomleft'
       , legend=c('representing trees, diverse set', 'represented tree')
       , col=c('blue','green')
       , pch=c(0,1)
       , cex = 1
       , bty='n')  
  return(list('div_rep_sf'=drsf, metric=attr(D,'metric'), 'lod'=lod))
}

#######################################################

# open first forest
load('data/nursery02/nursery02_01.rda')

DM <- doc[[1]]$DM

# run for both dissimilarities, produces figures 6.8

metric <- 'd0'
#metric <- 'd1'
{
  D <- DM[[metric]]
  attributes(D)
  attr(D,'metric') <- metric
}

m <- cmdscale(D, eig = TRUE, k = 2)

metric
attr(D,'metric')

par(mar=c(4,4,1,1)+0.1)

#plot_mds_forest(D,m)
stress(D,m)

plot_mds_representation(D,m,0.1)
mds_div_rep_sf <- plot_mds_div_rep_sf(D,m,0.1)

sf <- mds_div_rep_sf$div_rep_sf
lod <- mds_div_rep_sf$lod ; lod

plot(x=as.dist(D[sf,sf]), y=dist(m$points[sf,]) , xlab='dissimilarity' , ylab='euclidean distance in plane')
X <- cbind(A=as.dist(D[sf,sf]),B=dist(m$points[sf,]) ) %>% data.frame
lmX <- lm(B~A, data=X)
abline(lmX$coefficients, col='grey')
lmX %>% summary

#### not in thesis
#### shows that the mds does not work with 500 trees - but works ok with 10 trees

# 500 trees are an awful lot to plot in 2 dimensions . Reduce to 10.
D_small <- D[1:10,1:10]
m_small <- cmdscale(D_small,k=2,eig=T)
# plot_mds_representation(D_small,m_small,0.1)
plot_mds_representation(D_small,m_small,0.1)
plot_mds_div_rep_sf(D_small,m_small,0.1)

stress(D_small,m_small)
plot(x=as.dist(D_small), y=dist(m_small$points))
lm(as.dist(D_small)~dist(m_small$points))$coefficients
# still a lot of stress, but the plot makes sense


