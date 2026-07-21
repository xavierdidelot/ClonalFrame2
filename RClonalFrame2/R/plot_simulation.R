#' Plotting simulated recombination events from cf2out objects
#' @param x Output from ClonalFrame2
#' @param showSubs Whether to show the individual substitutions on each branch
#' @param ... Additional parameters are passed on
#' @return Plot of ClonalFrame2 simulation
plot_simulation = function(x,showSubs=F,...) {
  newick=x$simulation_parameters$tree
  tree=plot_tree_aligned(newick,x$simulation_parameters$length)

  #Show substitutions
  if (showSubs) {
    ancseqs=x$ancestral_sequences
    for (i in (1:Ntip(tree))) {
      father=which(tree$edge[,2]==i)
      if (tree$edge.length[father]==0) father=which(tree$edge[,2]==tree$edge[father,1])
      father=tree$edge[father,1]
      father2=which(tree$edge[,1]==father & tree$edge.length==0)
      w1=which(ancseqs$patmatrix$branch==tree$tip.label[i])
      if (length(father2)==0) w2=which(ancseqs$patmatrix$branch==sprintf('NODE%d',Ntip(tree)))
      else {
        father2=tree$edge[father2,2]
        w2=which(ancseqs$patmatrix$branch==tree$tip.label[father2])
      }
      s1=strsplit(ancseqs$patmatrix[w1,2],'')[[1]]
      s2=strsplit(ancseqs$patmatrix[w2,2],'')[[1]]
      s1=s1[ancseqs$patindex+1]
      s2=s2[ancseqs$patindex+1]
      diffs=which(s1!=s2)
      points(diffs,rep(i,length(diffs)),col='black',pch=4)
    }
  }

  #Show recombination events
  recevents=x$simulation_parameters$recombination_events
  for (i in (1:Ntip(tree))) {
    w=which(recevents$branch==tree$tip.label[i])
    if (length(w)==0) next
    segments(recevents[w,2],i,recevents[w,3],i,col='blue',lwd=2)
  }

}
