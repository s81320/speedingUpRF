I first ran code/nursery/01.R with creating the d1 dissimilarity matrix
with the full training data as argument.
a different 2nd argument like this:
DM$d1 <- createDMd1(forest=rg$forest, dft=data.train)

with data.train the Cleveland simulation, the training data for the forest.
This created 20 files in data/nursery.

BUT: This did contain doubles, observations repetedly sampled to be in the simulation. 
And doubles were always placed in the same terminal node. But that was not, what we
were interested in: to count the doubles in the simulation.

So I removed doubles and updated the d1 dissimilarity matrix as

DM$d1 <- createDMd1(forest=rg$forest, dft=Cleve[cr[[i]],])

with cr the (re)sample, so dft now was the set of unique observations, a subset
of the original Cleveland data set. This was done with the script update_d1.R and it
created the 20 files with updated d1 dissimilarity matrices in data/nursery02.

I also corrected a typo: In the first run I saved the training data under a misspelled name

'bootstapped training data'=data.train (r missing in bootstRapped)

code/nursery/01.R then was corrected.

When now running code/nursery/01.R it will create the d1 dissimilarity in the updated way
and the script update_d1.R is no longer needed.

################################################################################

Script sim1_with_inbag_info.R created data/sim1_with_inbag_info.rda which contains
basically the same simulations, forests and dissimilarity matrices as 
data/nursery02/nursery02_01.rda with only one change: The ranger forests were
created with the option keep.inbag=TRUE. This was needed to access the exact 
training data per tree: the original Cleveland observations that went into the 
simulation and then into training the tree.