# exploratory data analysis

# chapter data set 
# chapter 1


{rm(list=ls())
load('data/data_SupMat4.rda')
df.list <- ls()
print(df.list)
}

# Cleveland for i=1
{i <- 1
print(paste('now working on',df.list[[i]]))
df <- get(df.list[[i]])[,1:11] # omit last column, target with different level names
attr(df,'data.set.name') <- df.list[[i]]
nums <- names(df)[Vectorize(function(i){is.numeric(df[,i])})(1:11)] # numerical features

summary(df) %>% print
df[,nums] %>% apply(2,(function(x) c(range(x),mean(x),sd(x))))
}

# visual check if numerical features are continuous or discrete
{par(mfrow = c(1, length(nums)))
  par(mar=c(4,2,1,1)+0.2)
  for(n in nums){
  hist(df[,n], main=n)
}
}

# plotting features as density (numerical) or barchart (categorical)
# density can be misleading: constant variable is plotted as nice unimodal density -> check for sd>0 before!
{filename <- paste('plots/intro/features-overview-',attributes(df)$data.set.name,'.png',sep='')
png(file=filename, width=800, height=600, unit='px')

#dev.new(width = 50, height = 12, unit = "px")
par(mfrow = c(3, 4))
par(mar=c(4,2,1,1)+0.2)
for(n in names(df)){
  if(is.numeric(df[,n])){
    
    plot(density(df[,n]), main='', xlab=n, ylab='', yaxt='n', xaxt='n', cex.lab=2.5)
#    hist(Cleve[,n], main='', xlab=n, ylab='', xaxt='n', cex.lab=2.5)
  }else{
    plot(df[,n], main='' , xlab=n, ylab='', xaxt='n', cex.lab=2.5)
  }
}
par(mfrow = c(1, 1))
dev.off()
}

plot(density(Cleve[,'Age']), main='', xlab=n, ylab='', xaxt='n', cex.lab=2.5)

"n <- 'RestingECG'
#n <- 'Chestpaintype'
Cl2 <- Cleve[,n]
levels(Cl2)
if(n=='Chestpaintype'){
  levels(Cl2) <-  c('Typ. ang.', 'Atyp. ang.','non-ang.','asympt.')
}

if(n=='RestingECG'){
  levels(Cl2) <-  c('Normal', 'ST-T w. a.','Left hyper.')
}

plot(Cl2, main=n)
"

# no perfect correlation, decent correlation, some close to uncorrelated

par(mfrow = c(1, 1))
pairs(df[,nums], cex.labels=2.5, xaxt='n', yaxt='n')
corrplot::corrplot(cor(df[,nums]), method='pie', type = 'upper', diag = F)

#### comparing data sets

{
  emptyplot <-  function(LABEL){
    # empty plot for row label
    plot(x=c(-1,1) , y=c(-1,1), type='n', xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
    text(0,0,LABEL, cex=2.5)
    }
  
  filename <- 'plots/intro/compare-features-over-datasets-2.png'
  png(file=filename, width=900, height=1400, unit='px')
  dfs <- df.list[1:3]
  #dev.new(width = 50, height = 12, unit = "px")
  par(mfrow = c(12, 4))
  par(mar=c(4,2,1,1)+0.2)
  
  emptyplot('')
  emptyplot('Cleveland')
  emptyplot('Hungary')
  emptyplot('Swiss')
  
  for(n in names(df[1:10])){
    
    # empty plot for row label
    emptyplot(n)
    
    if(n %in% nums){
      # for numerical features: create a common xlim for all data sets
      lapply(dfs, function(x) get(x) %>% .[,n]) %>% 
        unlist %>% 
        range -> xlim
    }#else{
    # for categorical features: create a common ylim for all data sets
    # one should also check if all cat features use the full range of levels.
    #lapply(dfs, function(x) table(get(x)[,n]) %>% max) %>% unlist %>% unname %>% max %>% c(0,.) -> ylim
    #}
    for(X in dfs){
      X %>% get %>% .[,n] -> X
      if(n %in% nums){
        hist(X, main='', xaxt='n', yaxt='n',xlab='', xlim=xlim)
      }else{
        plot(X, main='' , xlab='', ylab='', xaxt='n', yaxt='n')
      }
      box()
    }
  }
  
  # empty plot for row label
  emptyplot('CAD')
  
  for(X in dfs){
    X <-  get(X)
    X <-  X[,'CAD']
    plot(X, main='' , xlab='', ylab='', xaxt='n', yaxt='n')
    box()

  }
  dev.off()
}
