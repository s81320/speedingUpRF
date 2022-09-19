# chapter methods
# plot a simple tree to visualize the nodeIDs in ranger
# plot its equivalent as a full binary tree

# chapter methods 
# section dissimilarity of trees
# figure 2.4 fig:intro:2trees:ranger and figure 2.5 fig:intro:2trees:fbt

source('code/source/plotTree.R')

load('data/nursery02/nursery02_01.rda')

plotTree1(rg= doc[[1]]$rg ,tri=1, plot.splitval = F)
plotTreeFb(rg= doc[[1]]$rg,1)

plot_el <- function(el, vlc=NULL, vl0=T){
  graph_from_edgelist(el) -> g1
  
  lo <- layout_as_tree(g1,root=1)
  
  if(is.null(vlc)){
    vertex.label.color <- 'blue'
  }else{
    vertex.label.color <- ifelse(vlc==0,'blue','grey')
    #vertex.label.color <- ifelse(vlc==0,'blue','green')
  }
  
  if(vl0){
    vertex.label <- 0:(max(el[,2])-1)
  }else{
    vertex.label <- 1:max(el[,2])
  }
  
  g1 %>% plot(layout=lo 
              , vertex.label.cex=2
              , vertex.label = vertex.label
              , vertex.label= 1:max(el[,2])
              , vertex.label.color=vertex.label.color
              , vertex.shape='none'
              , edge.color ='red'
              , edge.arrow.size=.4
              , rescale = T
              , asp = 0
  )
  
  
}

# 1st tree, ranger node.IDs , indicating terminal nodes
{
  el <- matrix(0,14,2)
el[1,] <- c(1,2)
el[2,] <- c(1,3)
el[3,] <- c(2,4)
el[4,] <- c(2,5)
el[5,] <- c(3,6)
el[6,] <- c(3,7)
el[7,] <- c(4,8)
el[8,] <- c(4,9)
el[9,] <- c(5,10)
el[10,] <- c(5,11) 
el[11,] <- c(6,12) 
el[12,] <- c(6,13)
el[13,] <- c(7,14)
el[14,] <- c(7,15)
vlc <- rep(0,14) # vlc : v_ertex l_abel c_olour
vlc[10:12] <- 1
plot_el(el,vlc)
}

# 2nd tree , ranger node.IDs , indicating terminal nodes
{
  el <- matrix(0,14,2)
  el[1,] <- c(1,2)
  el[2,] <- c(1,3)
  el[3,] <- c(2,4)
  el[4,] <- c(2,5)
  el[5,] <- c(3,6)
  el[6,] <- c(3,7)
  el[7,] <- c(4,8)
  el[8,] <- c(4,9)
  el[9,] <- c(5,10)
  el[10,] <- c(5,11) 
  el[11,] <- c(6,12) 
  el[12,] <- c(6,13) 
  el[13,] <- c(8,14) 
  el[14,] <- c(8,15) 
  vlc <- rep(0,14)
  vlc[c(7,10,11,15)] <- 1
  plot_el(el,vlc)
}


# tree 1 as full binary tree (height 5)
{
  hi <- 5
  el <- matrix(0,2^hi-2,2)
  ct <- 1
  for(n in 1:(2^(hi-1)-1)){
    el[ct,] <- c(n,2*n)
    el[ct+1,] <- c(n,2*n+1)
    ct <- ct+2
  }
  vlc <- ifelse(plotTreeFb(rg= doc[[1]]$rg,1)$sb.info[1:(2^hi-1)]==0,1,0)
  plot_el(el, vlc, vl0=F)
  plotTreeFb(rg= doc[[1]]$rg,tri=1)$sb.info[1:(2^hi-1)]
  }


# tree 2 as full binary tree (height 5)
{
  hi <- 5
  el <- matrix(0,2^hi-2,2)
  ct <- 1
  for(n in 1:(2^(hi-1)-1)){
    el[ct,] <- c(n,2*n)
    el[ct+1,] <- c(n,2*n+1)
    ct <- ct+2
  }
  vlc <- ifelse(plotTreeFb(rg= doc[[1]]$rg,2)$sb.info[1:(2^hi-1)]==0,1,0)
  plot_el(el, vlc, vl0=F)
  plotTreeFb(rg= doc[[1]]$rg,tri=2)$sb.info[1:(2^hi-1)]
}

