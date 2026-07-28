# Analysis of many simulated datasets
library(RClonalFrame2)
rm(list=ls())
n=100#Number of sequences
l=1e6#Genome length
burnin=0.5
allres=c()
nams=c('rho','delta','nu')
for (ind in 1:100) {
  filename=sprintf('long%d.RData',ind)
  if (!file.exists(filename)) next
  load(filename)
  res=c()
  for (nam in nams) {
    #true value
    res=c(res,get(nam))
    #inference without recombination
    iters=max(1,floor((nrow(x0$iterations)-1)*burnin)+2):nrow(x0$iterations)
    v=x0$iterations[iters,nam]
    v=sort(v)
    vals=c(mean(v),v[pmax(1,floor(length(v)*c(0.025,0.975)))])
    res=c(res,vals)
    #inference with recombination
    iters=max(1,floor((nrow(x$iterations)-1)*burnin)+2):nrow(x$iterations)
    v=x$iterations[iters,nam]
    v=sort(v)
    vals=c(mean(v),v[pmax(1,floor(length(v)*c(0.025,0.975)))])
    res=c(res,vals)
  }
  allres=rbind(allres,res)
}
pdf('ana_simu_multi.pdf',12,5)
ranges=matrix(c(-2,0,2,3,-2,-1),3,2,byrow = T)
par(mfrow=c(1,3))
for (para in 1:3) {
  plot(-1,-1,xlim=ranges[para,],ylim=ranges[para,],xlab='',ylab='',pch=20)
  points(log10(allres[,(para-1)*7+1]),log10(allres[,(para-1)*7+5]))
  lines(ranges[para,],ranges[para,])
  for (i in 1:nrow(allres))
    lines(c(log10(allres[i,(para-1)*7+1]),log10(allres[i,(para-1)*7+1])),c(log10(allres[i,(para-1)*7+6]),log10(allres[i,(para-1)*7+7])))
}
dev.off()
browseURL('ana_simu_multi.pdf')
