#' Testing SPR in C++
#' @param n Number of sequences
#' @param niter Number of iterations
#' @returns Matrix of results
#' @export
test_spr=function(n=4,niter=1000)
{
  system(sprintf('../build-debug/tests/test_spr %d %d > res.txt',n,niter))
  trees=read.tree('res.txt')
  system('rm res.txt')
  clusters=rep(0,length(trees))
  reps=list()
  k=0
  for (i in seq_along(trees)) {
    assigned=FALSE
    for (j in seq_along(reps))
      if (RF.dist(trees[[i]], reps[[j]]) == 0) {clusters[i]=j;assigned=TRUE;break}
    if (!assigned) {k=k+1;reps[[k]]=trees[[i]];clusters[i]=k}
  }
  counts=table(clusters)
  print(counts)
  print(chisq.test(counts))
  allres=trees
  return(allres)
}
