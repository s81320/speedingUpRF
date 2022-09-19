# generates ranger forest from 1. simulation with parameter keep.inbag=T 
# is the same forest as in nursery / nursery02, with additional info (inbag observations per tree!)

# saves the forest 
# it is needed in the introductory chapter methods to present performance of individual trees
# and in review of methods to compare Chipman 1, 2 and Meiner forests


# file data/sim1_with_inbag_info.rda is needed in scripts ensemble-magic.R , compare-models.R
# this rda-file is created with this script
# if it already exists, it is not overwritten, but saved as ...*.rda


load('data/data_SupMat4.rda') # loads the data sets Cleve, Hung, Swiss, VA

# required for create_prob()
source('code/source/sim-prep-v2.R')
  
{
# block copied from code/nursery/01.R
# may not be changed! Creates exactly the same (first 3) trees as in first simulation
# but with added keep.inbag=T

Cleve_enriched <- create_prob(Cleve)
nBs <- 50
  
seed <- 13
set.seed(seed)
cr<-createResample(Cleve_enriched$CAD
                     , times = nBs
                     , list = T)
  
ct <- 1
i <- 1
data.train <- new_bootstrap(Cleve_enriched , cr[[i]])[,-12] # no probabilities
  
rg <- ranger(CAD~.
               , data = data.train 
               , num.trees = 500
               , replace = F 
               , mtry= 3 
               , importance = 'impurity'
               , probability = T 
               , min.node.size = 13 
               , keep.inbag = T # NEW # only difference to code/nursery/01.R
)
  
info <-  list('created with script file'='code/nursery/sim1_with_inbag_info.R' 
                , 'seed'=seed
                , 'nBs'=nBs 
                , 'date'=Sys.time() 
                )
  
doc2 <- list( 'resample'=cr[[i]]
                , 'bootstrapped training data'=data.train
                , 'rg'=rg
                )
}

filename <-  'data/sim1_with_inbag_info.rda'  
if(filename %in% list.files(recursive = T)){
  print(paste(filename , 'already exists. File is saved with an * added to the name.'))
  filename <- strsplit(filename,split=".rda") %>% paste('*.rda', sep="")
  save(info, doc2, file=filename)
}else{
  save(info, doc2, file=filename)
}


