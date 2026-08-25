#' Comparison of likelihood functions
#' @param n Number of sequences
#' @param l Length of sequences
#' @param mu Mean of branch lengths
#' @param niter Number of iterations
#' @param eps Standard deviation of normal noise to apply to distances
#' @param showProgress Whether to show a progress bar
#' @param showPlot Whether to show a plot of the results
#' @param showTimes Whether to show the time spent by each method
#' @returns Matrix of results
#'
comparison=function(n=20,l=10000,mu=0.01,niter=100,eps=0.01,showProgress=T,showPlot=T,showTimes=T) {
  tree=rtree(n,br=function (n) rexp(n,rate=1/mu))
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
  if (showPlot && !requireNamespace("GGally", quietly = TRUE)) message("Package 'GGally' is required to use showPlot")
  if (showPlot &&  requireNamespace("GGally", quietly = TRUE)) {p=GGally::ggpairs(allres);print(p)}
  return(allres)
}


#' Adding random noise to a phylogenetic tree
#' @param tree Tree to randomise
#' @param eps Standard deviation of normal noise to apply to distances
#' @returns Randomised tree
#'
tree_noise=function(tree,eps=0.01) {
  d=cophenetic.phylo(tree)
  d=as.dist(d)
  d=abs(d+rnorm(length(d),0,eps))
  p=nj(d)
  p$edge.length=pmax(p$edge.length,1e-6)
  p=phangorn::midpoint(p)
  return(p)
}

#' Comparison of likelihood functions in C++
#' @param n Number of sequences
#' @param l Length of sequences
#' @param mu Mean of branch lengths
#' @param niter Number of iterations
#' @param showPlot Whether to show a plot of the results
#' @returns Matrix of results
#'
comparisonCpp=function(n=20,l=10000,mu=0.01,niter=100,showPlot=T) {
  restxt=tempfile()
  system(sprintf('../build-release/comparison %d %d %f %d %d > %s',n,l,mu,niter,round(runif(1,min=0,max=1e9)),restxt))
  allres=read.table(restxt,sep='\t',fill=T)
  allres=allres[seq(2,nrow(allres),2),]
  allres[,1]=as.numeric(allres[,1])
  methods=c('pruning','joint','parsimony')
  colnames(allres)<-methods
  if (showPlot && !requireNamespace("GGally", quietly = TRUE)) message("Package 'GGally' is required to use showPlot")
  if (showPlot &&  requireNamespace("GGally", quietly = TRUE)) {p=GGally::ggpairs(allres);print(p)}
  return(allres)
}
