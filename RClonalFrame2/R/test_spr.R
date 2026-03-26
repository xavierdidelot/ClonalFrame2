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
  allres=trees
  return(allres)
}
