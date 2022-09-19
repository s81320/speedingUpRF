# 2.9.2022

# compare test result
# which model is best?
# plot success (=logloss) over size (=number of trees)

rm(list=ls())

library(dplyr)

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

# for tested Chipman based model
{filename <- 'data/testNewMethods/test-all-methods-20files.rda'
        load(filename)
        r1 <- res$res
        
        r1 %>% 
                group_by(method, purpose, dissimilarity, p.representation) %>% 
                summarize(m.size=mean(size)
                        #, sd.size=sd(size)
                        , m.alpha=mean(alpha)
                        , m.chf=mean(logloss.chf)
                        #, sd.chf=sd(logloss.chf)
                        , m.bc=mean(logloss.base.case)
                ) -> tb2
}


results_plot <-  function(tb1, tb3){
#tb1 <- tbA ; tb3 <- tbB
        ylim <- range(rbind(tb1[,c('m.chf','m.bc')],tb3[,c('m.chf','m.bc')]))
        
        col.map <- c('chip1'='blue'
                     ,'chip2'='orange'
                     ,'chip1_simplified'='green'
                     ,'meiner'='red')
        legend.map <- c('chip1'='Chipman 1'
                        ,'chip2'='Chipman 2'
                        ,'chip1_simplified'='simplified Chipman 1'
                        ,'meiner'='Meiner')
        dissim.map <- c('d0'=1 , 'd1'=2, 'd2'=3, 'sb'=4)
        
        plot(tb1[,c('size','m.bc')],  type='l', col='grey', ylim=ylim, ylab='logloss')
        points(tb1[,c('size','m.chf')],  type='l' , lty=2 , col='grey')
        
        points(tb3[,c('m.size','m.chf')]
               , col=col.map[tb3$method]
               , pch=as.numeric(factor(tb3$dissimilarity, levels=names(dissim.map)))
        )
        dssms <- levels(as.factor(tb3$dissimilarity)) # dissimilarities
        mthds <- levels(as.factor(tb3$method)) # methods
        legend( 'topright'
        #         'bottomleft'
               , legend = c(legend.map[mthds],dssms)
               , col=c(col.map[mthds],rep(1,length(dssms)))
               , pch=c(rep(15,length(mthds)),dissim.map[dssms])
               , bty='n'
               #, cex=0.9
        )
        legend('bottomleft'
               , legend=c('regular small forest','high performers')
               , col='grey'
               , lty=c(1,2)
               , bty='n'
               )
}


# full view
names(tb2)
results_plot(tb1,tb2)

# zoom in on 5 to 10
{
        tbA <- tb1[tb1$size %in% 5:10,]
mask <- (tb2$m.size >=5) & (tb2$m.size <=10)
tbB <- tb2[mask,c('m.size','m.chf','m.bc','method','dissimilarity')]
results_plot(tbA, tbB)
}

# zoom in on 40 to 50
{tbA <- tb1[tb1$size %in% 40:50,]
mask <- (tb2$m.size >=40) & (tb2$m.size <=50) & (tb2$purpose=='mit50')
tbB <- tb2[mask,c('m.size','m.chf','m.bc','method','dissimilarity')]
results_plot(tbA, tbB)
}
