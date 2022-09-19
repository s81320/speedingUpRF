# test high performers against regular small forests

rm(list=ls())

library(dplyr)

# compare two sets of high performers and base cases
{filename <- 'data/testNewMethods/test-hp-5to50-20files.rda'
  load(filename)
  r1 <- res$res
  
  r1 %>% 
    group_by(method, size) %>% 
    summarize( m.chf=mean(logloss.chf)
               , sd.chf=sd(logloss.chf)
               , m.bc=mean(logloss.base.case)
               , m.size=mean(size)
               , m.alpha=mean(alpha)
    ) -> tb1
  tb1$overshoot <- tb1$m.chf - tb1$m.bc
  
  filename <- 'data/testNewMethods/test-hp-5to50-20files-2.rda' # selected trees starting at 350+1 for the basecase
  load(filename)
  r2 <- res$res
  
  r2 %>% 
    group_by(method, size) %>% 
    summarize( m.chf=mean(logloss.chf)
               , sd.chf=sd(logloss.chf)
               , m.bc=mean(logloss.base.case)
               , m.size=mean(size)
               , m.alpha=mean(alpha)
    ) -> tb2
  tb2$overshoot <- tb2$m.chf - tb2$m.bc
  
  # overshoot evolves the same way for tb1 and tb2
  plot(tb1$overshoot, type='l')
  points(tb2$overshoot, type='l', lty=2)
  
  # base case and high performers, tb1
  plot(tb1$m.chf, type='l')
  points(tb1$m.bc, type='l', lty=2)
  #zoom to 40 to 50 : basecase has lower logloss than high performers
  plot(tb1[37:46,]$m.chf, type='l')
  points(tb1[37:46,]$m.bc, type='l', lty=2)
  
  # base case and high performers, tb2
  plot(tb2$m.chf, type='l')
  points(tb2$m.bc, type='l', lty=2)
  #zoom to 40 to 50 : basecase has lower logloss than high performers
  plot(tb2[37:46,]$m.chf, type='l')
  points(tb2[37:46,]$m.bc, type='l', lty=2)

}

# for high performers
{filename <- 'data/testNewMethods/test-hp-5to50-20files.rda'
  load(filename)
  r1 <- res$res
  r1$size <- r1$size %>% 
    as.character %>%
    as.factor # works for the parameter p.target.size , gets hickups if target size is not met
  #r1[r1$purpose=='mit5to10',] %>% 
  r1 %>% 
    group_by(method, size) %>% 
    summarize( m.chf=mean(logloss.chf)
               #, sd.chf=sd(logloss.chf)
               , m.bc=mean(logloss.base.case)
               #, m.size=mean(size)
               #, sd.size=sd(size)
               # m.alpha=mean(alpha)
    ) -> tb1
  
  tb1$size %>% 
    as.character %>% # changes to levels first
    as.numeric -> tb1$size
  
  tb1 <-  tb1[order(tb1$size),]
}

# evolution of overshoot as sizes grow
{filename <- 'data/testNewMethods/test-hp-5to500-20files.rda'
  load(filename)
  r1 <- res$res
  X1 <- r1$logloss.chf - r1$logloss.base.case
#for(n in 5:50){
for(n in unique(r1$size)){  
r2 <- r1[r1$size==as.character(n),c('logloss.chf','logloss.base.case')]
X2 <- r2$logloss.chf - r2$logloss.base.case
hist(X2
     , breaks = seq(min(r1$logloss.chf - r1$logloss.base.case), 0.05+max(r1$logloss.chf - r1$logloss.base.case),0.05)
     , main=paste('size:',n,', mean overshoot: ',round(mean(X),6)) 
     , xlab='logloss overshoot'
     , xlim=range(X1)
     , ylim=c(0,1000)
     )
}
}
