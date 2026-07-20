#' Extract a given iteration from cf2out objects
#' @param x Output from ClonalFrame2
#' @param iteration Index of the iteration to extract (default is the last iteration)
#' @return Object of class cf2out corresponding to the extracted iteration
extract_iteration = function(x,iteration) {
  stopifnot(inherits(x, "cf2out"))
  if (missing(iteration)) iteration=nrow(x$iterations)

  #Extract properties of iteration
  rho=x$iterations[iteration,'rho']
  delta=x$iterations[iteration,'delta']
  nu=x$iterations[iteration,'nu']
  tree=read.tree(text=x$iterations[iteration,'tree'])
  likelihood=x$iterations[iteration,'likelihood']

  #Recover data
  ancseqs=x$ancestral_sequences
  a=ancseqs$patmatrix$data[1:Ntip(tree)]
  names(a)=ancseqs$patmatrix$branch[1:Ntip(tree)]
  l=as.list(a)
  seqs_split=lapply(l, function(x) strsplit(x, "")[[1]])
  seqs=phyDat(seqs_split, type = "DNA")
  attr(seqs,'index')=ancseqs$patindex+1
  attr(seqs,'weight')=tabulate(attr(seqs,'index'), nbins = attr(seqs,'nr'))

  #Perform ancestral state reconstruction
  pml=phangorn::pml(tree,seqs,model='JC')
  anc=phangorn:::joint_pml(pml)

  #Compute likelihood and check it is the same as in C++
  #d=likelihood_anc_rec(tree,seqs,anc,rho=rho,delta=delta,nu=nu)
  #print(c(likelihood,d))

  #Inference of recombined regions
  x2=decode_recombination(tree,seqs,anc,rho,delta,nu)
  print(c(likelihood,x2))
  return(x2)
}

#' Plotting tree and recombination events of a given iteration from cf2out objects
#' @param x Output from ClonalFrame2
#' @param iteration Index of the iteration to plot (default is the last iteration)
#' @param ... Additional parameters are passed on
#' @return Plot of ClonalFrame2 result
plot_iteration = function(x,iteration,showSubs=F,...) {
  stopifnot(inherits(x, "cf2out"))
  if (missing(iteration)) iteration=nrow(x$iterations)
  x2=extract_iteration(x,iteration)
  plot_recombination(x2,showSubs,...)
}

