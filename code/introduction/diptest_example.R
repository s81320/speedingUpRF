# checking out the Hartigans' dip test
# 5.5.2022
# example for the diptest applied to higher dimensional data
# in thesis : figure fig:unimod:1 , figure 2.7
# code / introduction / diptest example.R


par(mar=c(4,4,3,1)+0.1)


#### create data as multimodal. 

#### the dip test is for 1 dimensional data
#### now: high dimensional data : 2 dim as example

# install.packages('MASS')
library(MASS)
mu1 <- c(-4,0)
mu2 <-  c(0,4)

sigma0 <- matrix(c(1,0,0,1),2,2)

set.seed(2)
data <- rbind(mvrnorm(10,mu1,sigma0),
              mvrnorm(10,mu2,sigma0))

par(mar=c(2,2,1,1)+0.1)

# plot to find prominent points that 'see' and do not 'see' clusters 
{
  plot(data, type='n' , xlab='' , ylab='')
  text(data, labels = 1:20, col=c(rep('green',10),rep('blue',10)))
  points(rbind(mu1,mu2),col='red',pch=0)
  legend('bottomright'
       , legend=c('mu1 , mu2','centered at mu1','centered at mu2') 
       , col=c('red','green','blue')
       , pch=c(0,1,1)
       , cex=1.2)
}
# RESULT : select 6 ('sees' clusters) and 10 (does not 'see' clusters)

# scatterplot in figure 3.7
{
  par(mar=c(2,2,1,1)+0.1)

  plot(data , xlab='' , ylab='' 
     , col=c(rep('green',10),rep('blue',10))
     , pch=1)
  points(rbind(mu1,mu2),col='red',pch=0)
  legend('bottomright'
       , legend=c(expression(mu[1]*','*mu[2]),expression('centered at '*mu[1]), expression('centered at '*mu[2]))
       , col=c('red','green','blue')
       , pch=c(0,1,1)
       , cex=1.2)
  text(x=data[c(6,10),1]+0.1 , y=data[c(6,10),2]+0.1 , labels = c('2','1'), col='black')
# applying the dip test when working with a dissimilarity matrix
}

# creating the dissimilarity matrix based on points (R^k)
dm <- data %>% dist %>% as.matrix
# dip test for all data points (each points distance to all other points)
doc.dip <-  rep(NA,nrow(dm))
for(i in 1:nrow(dm)){
  doc.dip[i] <- dip.test(dm[i,-i])$p.value 
} 

################################################################################

# which data points 'see' clusters?
which(doc.dip<0.05)
# doc.dip[order(doc.dip)]
# order(doc.dip) 
# plot(doc.dip[order(doc.dip)], type='l') 
{
  pch <-  ifelse(doc.dip<0.05,16,1)
  plot(data , xlab='' , ylab='' 
     , col=c(rep('green',10),rep('blue',10))
     , pch=pch)
  points(rbind(mu1,mu2),col='red',pch=0)
  legend('bottomright'
       , legend=c(expression(mu[1]*','*mu[2]),expression('centered at '*mu[1]), expression('centered at '*mu[2]), 'sees clusters')
       , col=c('red','green','blue', 'black')
       , pch=c(0,1,1,16)
       , cex=1.2)
}

################################################################################

# histograms for bottom row in figure 3.7
hist(dm[6,-6],breaks=seq(0,8,0.5), main='', xlim=c(0,8), yaxt='n')
axis(side=2, at=seq(0,4))
plot(density(dm[6,-6]),main='')
dip.test(dm[6,-6])

hist(dm[10,-10],breaks=seq(0,8,0.5), main='', xlim=c(0,8) , yaxt='n', ylim=c(0,4))
axis(side=2, at=seq(0,4))
plot(density(dm[10,-10]),main='')
dip.test(dm[10,-10])

# discrete variables are always multimodal !!
# we can get a density - but a density is not for discrete variables. It does not make sense!!

