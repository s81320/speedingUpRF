# testing the Shannon-Banks dissimilarity

source('code/source/distance-matrices.R')

# defining full binary trees tri1 and tri2
# and counting the differences with function sb
# function sb is essential in createDMsb

# compare results with difference in hand drawn trees
# test case : 2 trees of same depth
tri1 <-  c(1
           ,2,3
           ,0,0,2,0)
# tri1 is 
#           1
#         /   \
#       2      3
#      / \    / \
#     0   0  2   0

tri2 <- c(2
          ,1,3
          ,2,0,0,0)
# tri2 is 
#           2
#         /   \
#       1      3
#      / \    / \
#     2   0  0   0

sb(tri1,tri2)
which(tri1!=tri2)
length(which(tri1!=tri2))
assertthat::assert_that(sb(tri1,tri2)==length(which(tri1!=tri2)))

# test case : 2 trees of different depth
# compare results with hand drawn trees
tri1 <-  c(1
           ,2,3
           ,0,0,2,0) # same as before 
tri2 <- c(2
          ,1,3
          ,2,0,0,0
          ,1,2,0,0,0,0,0,0)
# tri2 is 
#             2
#          /     \
#        1          3
#      /  \       /   \
#     2    0     0     0
#    / \  / \   / \   / \
#   1  2 0  0  0   0 0   0
assertthat::assert_that(sb(tri1,tri2)==sb(tri2,tri1)) # symmetry
assertthat::assert_that(sb(tri1,tri2)==6)

# test case : 2 trees of different depth, more than 1 layer difference
# compare results with hand drawn trees
tri1 <-  c(1,2,3,0,0,2,0) # same as before
tri2 <- c(2 # 1 
          ,1,3 # 2
          ,2,rep(0,3) # 4
          ,1,2,rep(0,6) # 8
          ,2,3,rep(0,14)) # 16
# tri2 is 
#                      2
#                  /       \
#                /           \ 
#              /               \
#             1                  3
#           /   \               /   \
#          /     \             /      \
#        2         0          0         0
#      /  \      /   \      /   \      /  \
#     1    2    0     0    0     0    0     0
#    / \  / \  / \   / \  / \   /  \ / \   /  \
#   2   3 0 0 0   0 0  0 0  0  0  0 0  0  0    0 
assertthat::assert_that(sb(tri1,tri2)==sb(tri2,tri1)) # symmetry
assertthat::assert_that(sb(tri1,tri2)==8) # 2 more non-zero entries in tri2 -> sb grows by 2, now is 8

################################################################################

source('code/source/plotTree.R')

# grow a small forest
rg <- ranger(Species~. , num.trees=3 , data=iris , max.depth = 5)

# test that treeInfo and the plotted trees agree
# maybe draw tree by hand from treeInfo to compare

treeInfo(rg,1)
plotTree1(rg, 1, main='tree 1', plot.splitval = F)
plotTreeFb(rg, 1, main='tree 1')

treeInfo(rg,2)
plotTree1(rg, 2, main='tree 2', plot.splitval = F)
plotTreeFb(rg, 2, main='tree 2')

treeInfo(rg,2) 
plotTree1(rg, 3, main='tree 3' , plot.splitval = F)
plotTreeFb(rg, 3, main='tree 3')


#########################################################

# the function rt2fbt (ranger tree to full binary tree) 
# is essential in createDMsb
# and also in the plot function plotTreeFb (plot the tree as a full binary tree).
# so if the tree is plotted correctly as a subset of a full binary tree
# then the correct transformation is made from a ranger tree to a full binary tree

# if trees are plotted correctly, 
# especially plotTreeFb, then rt2fbt is successfully tested

################################################################################

createDMsb(rg$forest)

# the following function is used in createDMsb
# full binary trees built from ranger trees
lapply(1:rg$forest$num.trees, function(i) rt2fbt(rg$forest,i) )

#################

# use return values of plotTree1 of plotTreeFb to calculate d0 by hand
createDMd0(rg$forest)
