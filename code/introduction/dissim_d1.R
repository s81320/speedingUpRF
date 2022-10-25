# needed in the introduction of d1 , example

load('data/nursery02/nursery02_01.rda')

rg <- doc[[1]]$rg
rs <- doc[[1]]$resample
data.train <- doc[[1]]$'bootstrapped training data'

P1 <- predict(rg,type = 'terminalNodes', data = data.train)$predictions[,1:3]

# number of terminal nodes in tree 1
tri <- 1
tN <- tNodes(rg$forest,tri) ; tN ; length(tN)
# fix the k-th terminal node
k <- 1
# observations in k-th terminal node for the tree tri (tree's number in the forest)
obsk <- which(P1[,tri]==tN[k]) ; obsk

pairs <- function(y) y*(y-1)/2
rs[obsk] %>% length ; unique(rs[obsk]) %>% length ; unique(rs[obsk]) %>% length %>% pairs

# where are these observations in another tree (tree 2)?
tri.other <- 2
tNodes(rg$forest, tri.other)
P1[unique(rs[obsk]),tri.other] 

# how many are (again) in the same terminal nodes in the other tree?
# terminal nodes ,for unique observations in terminal node k in tree 1, in tree 2
P1[unique(rs[obsk]),tri.other]
P1[unique(rs[obsk]),tri.other]  %>% table 
# number of pairs in terminal nodes in other tree
P1[unique(rs[obsk]),tri.other]  %>% table  %>% pairs
# total number of pairs
P1[unique(rs[obsk]),tri.other]  %>% table  %>% pairs %>% sum
# 2 pairs are (again) in the same terminal nodes

# rate at which pairs in one terminal node of a selected tree are again in the same terminal nodes in another tree
(P1[unique(rs[obsk]),tri.other]  %>% table  %>% pairs %>% sum) / (unique(rs[obsk]) %>% length %>% pairs)
# is a way to be similar, in this case: little similarity
# measure it as being dissimilar
1-(P1[unique(rs[obsk]),tri.other]  %>% table  %>% pairs %>% sum) / (unique(rs[obsk]) %>% length %>% pairs)
# VVVEEERRRYYY dissimilar

# check: if we remove uniqueness, the doubles should also be in the same 
# terminal nodes in the other tree , the total number of pairs increases
P1[rs[obsk],tri.other]  %>% table  %>% pairs %>% sum
