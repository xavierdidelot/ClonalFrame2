#' Plotting simulated recombination events from cf2out objects
#' @param x Output from ClonalFrame2
#' @param showSubs Whether to show the individual substitutions on each branch
#' @param ... Additional parameters are passed on
#' @return Plot of ClonalFrame2 simulation
plot_simulation = function(x,showSubs=F,...) {
  par(xpd=NA,mai=c(1,0.2,1,0.2),las=1)
  layout(matrix(c(1,2,2,2,2),nrow=1))
  genome_length=x$simulation_parameters$length
  newick=x$simulation_parameters$tree
  recevents=x$simulation_parameters$recombination_events
  ancseqs=x$ancestral_sequences

  #Add leaves corresponding to nodes
  repeat {
    new_newick = gsub("\\(([^(),]+),([^(),]+)\\)([^,:()]+):", "{\\1@\\3:0@\\2}:", newick)
    if (new_newick == newick) break
    newick = new_newick
  }
  newick = gsub("\\{","(",newick)
  newick = gsub("\\}",")",newick)
  newick = gsub("@",",",newick)
  tree=read.tree(text=newick)
  tree$node.label=NULL

  x=plot(tree,align.tip.label = T,tip.color='white',node.color='white',edge.color='white')
  #Drawing internal edges manually, otherwise they are not aligned with extra leaves
  disttoroot=dist.nodes(tree)[Ntip(tree)+1,]
  for (i in (1:Ntip(tree))) {
    father=which(tree$edge[,2]==i)
    if (tree$edge.length[father]==0) father=which(tree$edge[,2]==tree$edge[father,1])
    father=tree$edge[father,1]
    start=disttoroot[father]
    end=disttoroot[i]
    segments(start,i,end,i)
    father2=which(tree$edge[,1]==father & tree$edge.length==0)
    if (length(father2)>0) segments(start,i,start,tree$edge[father2,2])
  }
  fromroot=tree$edge[which(tree$edge[,1]==(Ntip(tree)+1)),2]
  for (i in 1:length(fromroot)) {
    w=which(tree$edge[,1]==fromroot[i] & tree$edge.length==0)
    if (length(w)>0)
    fromroot[i]=tree$edge[w,2]
  }
  segments(0,min(fromroot),0,max(fromroot))
  axisPhylo(1,backward = F)
  plot(c(),c(),xaxs="i",xlim=c(0,genome_length),ylim=c(1,Ntip(tree)),xlab='Genomic Position',ylab='',yaxt='n')
  labs=tree$tip.label
  labs=gsub("NODE[0-9]+","",labs)
  axis(2,at=1:Ntip(tree),labels=labs,gap.axis=-100)

  #Show substitutions
  if (showSubs) {
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
  for (i in (1:Ntip(tree))) {
    w=which(recevents$branch==tree$tip.label[i])
    if (length(w)==0) next
    segments(recevents[w,2],i,recevents[w,3],i,col='blue',lwd=2)
  }

}
