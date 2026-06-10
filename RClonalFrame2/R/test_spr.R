#' Testing SPR in C++
#' @param n Number of sequences
#' @param niter Number of iterations
#' @returns List of sampled trees
#' @export
test_spr=function(n=4,niter=1000)
{
  seed=round(runif(1,min=0,max=1e9))
  system(sprintf('../build-debug/test_spr %d %d %d > res.txt',n,niter,seed))
  trees=read.tree('res.txt')
  system('rm res.txt')
  clusters=rep(0,length(trees))
  reps=vector('list',length(trees))
  k=0
  for (i in seq_along(trees)) {
    assigned=FALSE
    if (k>0) for (j in 1:k)
      if (phangorn::RF.dist(trees[[i]], reps[[j]]) == 0) {clusters[i]=j;assigned=TRUE;break}
    if (!assigned) {k=k+1;reps[[k]]=trees[[i]];clusters[i]=k}
  }
  counts=table(clusters)
  print(counts)
  print(chisq.test(counts))
  invisible(trees)
}
