library(RClonalFrame2,quietly=T)
library(tictoc)
data(simu)
set.seed(0)
method='joint'
niter=1000

if (F) {
  tic()
  r<-experimental_mcmc(seqs,method=method,niter=niter)
  toc()
  print(r$lik_trace[length(r$lik_trace)])#to check behaviour change
  plot(r$lik_trace,type='l')
} else {
  set.seed(0)
  p=profvis::profvis(r<-experimental_mcmc(seqs,method=method,niter=niter))
  htmlwidgets::saveWidget(p, "profile.html")
  browseURL("profile.html")
}
