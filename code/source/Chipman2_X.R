# another Chipman forest: simplification of Chipman 1 by clustering regardless!

grow_chipForest_1_simplified <-  function(dm, oLL, parameter, output=F){
  #' grows the Chipman 1 forest ('core technique' in the paper by Chipman, George, McCollouch) enforcing clustering with a fixed number of clusters
  #'
  #' needs data.test in parent environment
  #' @param dm full dissimilarity matrix to calculate the level of representation via quantils 
  #' @param parameter needs to have cutoff, selection, sizeSF
  #' @param output if TRUE the function outputs also ...
  
  #print(dim(dm))
  #print(pa)
  
  sizeSF <-  parameter$sizeSF
  # level of representation , remove the diagonal, it has only 0 values
  alpha <- quantile(dm[upper.tri(dm)],parameter$cutoff) # diagonal with all values 0 not included
  #oLL <-  order(LL)
  
  represented <-  function(I){
    #' testing for representation of best I trees (R= dense representing forest = {1,..,I}) at level alpha
    #' 
    #' uses N, dm2 and alpha from parent environment. i.e. function calc_chipForest_1
    #' for I >1 , 
    #' representation by one tree alone can happen only at level of largest dissimilarity
    #' no representation by almost all trees (I==N-1) only if last tree, 
    #' index N, is far away from all others, further than cutoff
    #assertthat::assert_that(class(I)=='numeric', 'argument I has to be numeric') %>% print
    #assertthat::assert_that( (I >0) , 'I must be at least 1 and at most one less than trees in forest') %>% print
    
    if(I==1){ # max
      dm2[1,2:N] %>% max %>%
        (function(x){x<=alpha}) 
    }else{
      if(I==N-1){ # min
        dm2[1:(N-1),N] %>% min %>%
          (function(x){x<=alpha})
      }else{ # min and max
        dm2[1:I,(I+1):N] %>%
          apply(2,min) %>% # not working for I = 1 or I=N-1
          max %>%
          (function(x){x<=alpha})
      }
    }
  }
  
  dm2 <-  dm[oLL,oLL]
  
  N <-nrow(dm2) # used to be min(nrow(dm2),N) , N would stop early as soon as N trees were selected into the dense representing subforest.
  # But this function is for the original approach of the Chipmen, and they would not have done early stopping.
  
  if(alpha<max(dm2)){
    I <- 2
    while((I<N) && !represented(I)){
      I <- I+1
    }
  }else{I <- 1}
  # exits with smallest I with represented(I) TRUE
  
  # corrected error 10.7.22: Used to be R <-  1:I
  R <-  oLL[1:I] # (dense) representing subforest , dense as all trees inbetween (from best oLL[1] to I-th best oLL[I]) are included in the subforest
 
  # I will NOT be changed further down , R may be reduced by clustering
  # attempt clustering only if the representing forest is large (larger than 10)
  if(length(R)>sizeSF){ # cluster only if you got more than you want!
    # attempting to cluster makes sense only for R relatively large. Larger than 10 trees, say. 
    # And number of cluster should be b/w 5 and 50? And not larger than I/4?
    # dissimilarity matrix needed only for R
    dm2 <-  dm[R,R] # using R # funny name in an R script ...
    # which 2 trees are meant in dm2[1,2] ?? the best and the second best , indexed oLL[1], oLL[2]
    
    # always cluster
    cluster::pam(x=dm2
        , k=sizeSF
        , diss=T
        , medoids = 1:sizeSF # start with the best trees as medoids
        , cluster.only=F # if T, are medoids returned? We need them!
        ) -> pam.clus
   
        # from each cluster select the tree with the smallest logloss (on OOB data)
    if(parameter$selection=='best'){
      lapply(1:sizeSF,function(i) which(pam.clus$clustering==i)[1] %>% oLL[.]) %>% # we take the first, because indices are ordered, first is smallest
        # oLL to go back to original index
        unlist -> cf # representing subset
        # if indices were not ordered, we'd do 
        #lapply(1:sizeSF,function(i) which(pam.clus$clustering==i) %>% min ) %>% unlist -> cf.alt
      }
        
    if(parameter$selection=='central'){
      pam.clus$medoids -> cf
    }
    
    }else{
      cf <-  R
    } # if more trees than needed for sizeSF
    
  
  doc.ret <- list('forest'=cf
                  , 'size.dense.representing.sf' = I
                  , 'alpha'= alpha
                  , 'call'=list(name='grow_chipForest_1b' , 'parameter'=parameter) )
  
  if(output){
    doc.ret$pam.sizes <- 1 # cluster sizes would be interesting to know
    doc.ret$quality <- 2 # some kind of quality measure would be good to know
  }
  
  return(doc.ret)
}

eval_chipForest_1_simplified <-  function(dm , forest, oLL, parameter, output=F){
  cf1b <-  grow_chipForest_1_simplified(dm , oLL, parameter)
  
  chipForest <-  cf1b$forest
  
  #print('eval_chipForest_1_simplified will return')
  #print(list(calcLogloss( subforest(forest, chipForest ), data.test)
  #           , cf1b$size.dense.representing.sf # I
  #           , length(chipForest)))
  return(list(calcLogloss( subforest(forest, chipForest ), data.test)
              , cf1b$size.dense.representing.sf # I
              , length(chipForest))) 
}
