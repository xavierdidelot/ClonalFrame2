#' Side-by-side plot without any event shown
#' @param newick Newick string of tree to show
#' @param genome_length Genome length
#' @returns Phylo object with extra leaves for nodes
#'
plot_tree_aligned=function(newick,genome_length,show.tip.label=T) {
  par(fig=c(0.05,0.3,0.1,0.9),mai=c(0,0,0,0), xaxs="i", yaxs="i",las=1,xpd=NA)

  #Add leaves corresponding to nodes
  original_newick=newick
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

  #Generate strings to identify leaves
  original_tree=read.tree(text=original_newick)
  strings=children_strings(original_tree)
  tree$tip.label_strings=tree$tip.label
  for (i in 1:length(tree$tip.label)) tree$tip.label_strings[i]=strings[tree$tip.label[i]]

  #Prepare right hand side of plot
  par(fig=c(0.23+0.07*show.tip.label,0.95,0.1,0.9),mai=c(0,0,0,0), xaxs="i", yaxs="i",las=1,xpd=NA,new=T)
  plot(c(),c(),xaxs="i",xlim=c(0,genome_length),ylim=c(1,Ntip(tree)),xlab='Genomic Position',ylab='',yaxt='n')
  labs=tree$tip.label
  labs=gsub("NODE[0-9]+","",labs)
  if (show.tip.label) axis(2,at=1:Ntip(tree),labels=labs,gap.axis=-100)
  else axis(2,at=1:Ntip(tree),labels=rep('',Ntip(tree)),gap.axis=-100)

  return(tree)
}
