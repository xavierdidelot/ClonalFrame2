#' Plotting tree and recombination events of a given iteration from cf2out objects
#' @param x Output from ClonalFrame2
#' @param iteration Index of the iteration to plot (default is the medoid)
#' @param gap Gap between sites in output (default is 100)
#' @param burnin Burnin proportion to apply if using medoid
#' @param show_truth Whether to show the truth for a simulated dataset
#' @param ... Additional parameters are passed on
#' @return Plot of ClonalFrame2 result
#'
plot_recombination = function(x,iteration,gap=100,burnin,show_truth=T,...) {
  stopifnot(inherits(x, "cf2out"))
  if (missing(burnin)) burnin=0.5
  if (missing(iteration)) {
    message("Computing medoid...")
    iteration=medoid(x,burnin)
  }
  message("Decoding for recombination...")
  r=decode_recombination(x,iteration,gap)

  message("Preparing plot...")
  newick=x$iterations$tree[iteration]
  genome_length=length(x$ancestral_sequences$patindex)
  tree_aligned=plot_tree_aligned(newick,genome_length)

  #Reorder rows of matrix r to match order of leaves in tree_aligned
  treeit=read.tree(text=newick)
  treeit_strings=children_strings(treeit)
  r2=matrix(0,nrow(r),ncol(r))
  for (i in 1:Ntip(tree_aligned)) {
    str=tree_aligned$tip.label_strings[i]
    w=which(treeit_strings==str)
    if (length(w)==0) next
    w=which(treeit$edge[,2]==w)
    r2[i,]=r[w,]
  }
  image(x=seq(1,genome_length,gap),y=1:nrow(r2),z=t(r2),useRaster=F,axes = F,add=T)

  if (show_truth && !is.null(x$simulation_parameters)) {
    #Show true recombination events
    truetree=read.tree(text=x$simulation_parameters$tree)
    perm=match(truetree$tip.label, treeit$tip.label)#reorder labels in the next 4 lines
    n_tips=Ntip(truetree)
    is_tip=truetree$edge[, 2] <= n_tips
    truetree$edge[is_tip, 2]=perm[truetree$edge[is_tip, 2]]
    truetree$tip.label=treeit$tip.label
    truetree_strings=children_strings(truetree)
    recevents=x$simulation_parameters$recombination_events
    for (i in (1:Ntip(tree_aligned))) {
      str=tree_aligned$tip.label_strings[i]
      truename=names(which(truetree_strings==str))
      w=which(recevents$branch==truename)
      if (length(w)==0) next
      segments(recevents[w,2],i,recevents[w,3],i,col='blue',lwd=2)
    }
  }
}

children_strings=function(tree) {
  n_tips=Ntip(tree)
  n_nodes=Ntip(tree) + Nnode(tree)
  desc=Descendants(tree, node = seq_len(n_nodes), type = "tips")
  mat=matrix(0L, nrow = n_nodes, ncol = n_tips,
                dimnames = list(seq_len(n_nodes), tree$tip.label))
  for (i in seq_len(n_nodes)) mat[i, desc[[i]]] <- 1L
  rownames(mat)=c(tree$tip.label,tree$node.label)
  strings=apply(mat, 1, paste, collapse = "")
  return(strings)
}
