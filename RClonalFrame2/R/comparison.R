#' Comparison of likelihood functions
#' @param n Number of sequences
#' @param l Length of sequences
#' @param mu Mutation rate per site per branch length
#' @param niter Number of iterations
#' @param eps Standard deviation of normal noise to apply to distances
#' @param showProgress Whether to show a progress bar
#' @param showPlot Whether to show a plot of the results
#' @param showTimes Whether to show the time spent by each method
#' @returns Matrix of results
#' @export
#'
comparison=function(n=20,l=10000,mu=0.01,niter=100,eps=0.01,showProgress=T,showPlot=T,showTimes=T) {
  tree=rtree(n)
  tree$edge.length=tree$edge.length*mu
  seqs=phangorn::simSeq(tree,l=l)
  methods=c('pruning','joint','parsimony')#can also add 'marginal'
  allres=matrix(NA,niter,length(methods))
  colnames(allres)<-methods
  workenv=new.env()
  if (showProgress) pb <- txtProgressBar(min=0,max=niter*length(methods),style = 3)
  tims=rep(NA,length(methods))
  names(tims)=methods
  trees=vector(mode='list',length=niter)
  for (j in 1:length(methods)) {
    gc(F)#call garbage collector to improve timing
    start=Sys.time()
    for (i in 1:niter) {
      if (showProgress) setTxtProgressBar(pb, i+(j-1)*niter)
      if (j==1) trees[[i]]=tree_noise(tree,eps)
        allres[i,j]=likelihood(trees[[i]],seqs,method=methods[j])
    }
    tims[j]=difftime(Sys.time(),start,units='secs')
  }
  if (showProgress) close(pb)
  if (showTimes) print(tims)
  if (showPlot) {p=GGally::ggpairs(allres);print(p)}
  return(allres)
}


#' Adding random noise to a phylogenetic tree
#' @param tree Tree to randomise
#' @param eps Standard deviation of normal noise to apply to distances
#' @returns Randomised tree
#' @export
#'
tree_noise=function(tree,eps=0.01) {
  d=cophenetic.phylo(tree)
  d=abs(d+rnorm(length(d),0,eps))
  p=nj(d)
  p$edge.length=pmax(p$edge.length,1e-6)
  return(p)
}

#' Comparison of likelihood functions in C++
#' @param n Number of sequences
#' @param l Length of sequences
#' @param mu Mutation rate per site per branch length
#' @param niter Number of iterations
#' @param eps Standard deviation of normal noise to apply to distances
#' @param showPlot Whether to show a plot of the results
#' @returns Matrix of results
#' @export
#'
comparisonCpp=function(n=20,l=10000,mu=0.01,niter=100,eps=0.01,showPlot=T) {
  system(sprintf('../build-debug/tests/comparison %d %d %f %d %f > res.txt',n,l,mu,niter,eps))
  allres=read.table('res.txt',sep=',')
  system('rm res.txt')
  methods=c('pruning','joint','parsimony')
  colnames(allres)<-methods
  if (showPlot) {p=GGally::ggpairs(allres);print(p)}
  return(allres)
}
