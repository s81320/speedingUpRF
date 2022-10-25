library(dplyr)

source('code/timbR.R')
source('code/source/distance-matrices.R')

load('data/data_SupMat4.rda')

load('data/nursery02/nursery02_01.rda') # loads DM with 4 dissimilarity matrices, one for each dissimilarity

rg <-  doc[[1]]$rg
DM <- doc[[1]]$DM
names(DM)
metric <- names(DM)[1] ; metric
dm <-  DM[[metric]]

Td0 <- measure_distances(rg, metric='splitting variables' )
(Td0 == DM[['d0']]) %>% all

microbenchmark::microbenchmark(measure_distances(rg, metric='splitting variables' )
                               , createDMd0(rg$forest)
                               , times=10
)
# timbR is 10 times faster, ups


# equality needs the test_data to be equal : only unique observations!
rs <- doc[[1]]$resample
Td1 <- measure_distances(rg, metric='terminal nodes', test_data = Cleve[unique(rs),1:11]) # start 11:50 bis 12:20 oder weniger
(Td1 == DM[['d1']]) %>% all

# some other test set : Hung
Td1 <-  measure_distances(rg, metric='terminal nodes' , test_data = Hung[1:100,])
Myd1 <-  createDMd1(rg$forest, dft=Hung[1:100,])
(Td1 == Myd1) %>% all

microbenchmark::microbenchmark(measure_distances(rg, metric='terminal nodes' , test_data = Hung[1:100,])
                               , createDMd1(rg$forest, dft=Hung[1:100,])
                               , times=10
)
# my implementation is 5 times faster (times=1)

# d2 in timbR only for classification trees
rg <-  ranger::ranger(formula='CAD~.' , data=Cleve[,1:11])
rg$forest$treetype
# maybe only for non pe forests?
Td2 <- measure_distances(rg, metric='prediction', test_data = Cleve[,1:11] ) # does not work for pe forests
Myd2 <- createDMd2(rg$forest,dft = Cleve[,1:11])
(Td2 == Myd2) %>% all

# make sure, that equality does not follow from some weired error, like both objects are NULL
Td2[1:2,1:10]
Myd2[1:2,1:10]

microbenchmark::microbenchmark('timbR measure distances'=measure_distances(rg, metric='prediction', test_data = Cleve[,1:11] )
                               , 'create DM d2'=createDMd2(rg$forest,dft = Cleve[,1:11])
                               , times=2
                               , check='identical' # supposed to check if both functions generate the same data
                               )


