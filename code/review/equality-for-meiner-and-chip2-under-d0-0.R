# compare Chipman 2 to Meiner forest 
# for dissimilarity d0 at representation parameter 0
# they should always be the same

# chapter review

load('data/nursery02/nursery02_01.rda') # loads doc

nReps <- length(doc)
  
for(i in 1:nReps){
    if(i%%10 ==0) print(paste(i/nReps , Sys.time()))
    
    DM <-  doc[[i]]$DM
    
    rg  <-  doc[[i]]$rg
    forest <-  rg$forest
    
    predict(forest
            , data=data.train
            , predict.all = T)$predictions[,2,] -> pp # dim 303 x 500
    
    # List of Loglosses : LL
    lapply(1:forest$num.trees
           , function(t){ 
             OOB <- which(rg$inbag.counts[[t]]==0) 
             pp[OOB,t] %>% 
               calcLogloss2( df=data.train[OOB,] ) %>% 
               unlist
           }
    ) %>% 
      unlist -> LL
    
    ch2 <- grow_chipForest_2(DM[['d0']], oLL=order(LL) , parameter=list(cutoff=0, sizeSF=500))
    mf <- grow_meinForest(DM[['d0']], LL=LL , parameter=list(cutoff=0, sizeSF=500))
    
    assertthat::assert_that(all(ch2$forest==mf$forest)) %>% print
    length(mf$forest) %>% print # to make sure they are not equal because both are of length 0 or length 500
}
