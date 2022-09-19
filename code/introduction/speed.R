# http://adv-r.had.co.nz/Profiling.html

# generated figure 3.2 labeled as fig:speed
# runtime for prediction over number of trees / size of forest

rm(list=ls())

# install.packages('microbenchmark')
library(microbenchmark)

library(dplyr)
library(ranger)

# my own code
source('code/source/subforest.R') # constructors for a sub-forest (class ranger.forest) and its hull (class ranger)

load('data/data_SupMat4.rda')

rg <- ranger(CAD~.
                    , data = Cleve[,1:11]
                    , num.trees = 500
                    , replace = F 
                    , mtry= 3 # default : 3
                    , importance = 'impurity'
                    , probability = T 
                    , min.node.size = 13 
)

# predictions with a ranger sub-forest of size i on the validation set
predSF<-function(i){predict(subforest(rg,1:i),data=Swiss)}

times<-100L # default
mb <- microbenchmark(predSF(5)
                     , predSF(50)
                     , predSF(100)
                     , predSF(150)
                     , predSF(200)
                     , predSF(250)
                     , predSF(300)
                     , predSF(350)
                     , predSF(400)
                     , predSF(450)
                     , predSF(500)
                     , times =times
)

p <- summary(mb)

View(p)
attr(p,'unit')

# par(mfrow=c(1,1))
xgrid<- c(5,seq(from=50, to=500 , by=50))

par(mar=c(4,4,1,1)+0.2)
plot(x=xgrid
     , y=p$mean
     #, main=paste('me(di)an execution times\n(times=',times,')')
     , ylab=paste('time in',attr(p,'unit'))
     , xlab='number of trees'
    # , ylim=c(6,18)
     , cex.axis=1.2
     , cex.lab=1.2
    , ylim=c(0,max(p$mean))
)

p$expr <- xgrid
lm1<-lm(mean~as.numeric(expr), data=p)
lm1 %>% summary
abline(coefficients(lm1), col='grey')

# summary(lm1)
lm1$coefficients

#### additional stuff ########################
##############################################

ggplot2::autoplot(mb)

par(mar=c(4,4,3,1)+0.2)
boxplot(mb
        , log = F
        , main=paste('outliers in execution times for ranger predictions\n(times=',times,')' ) # many outliers!
)

