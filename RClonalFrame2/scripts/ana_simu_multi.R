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

if (is.null(allres)) allres=matrix(1,100,100)
pdf('ana_simu_multi.pdf',10,4)
ranges=matrix(c(-2,0,2,3,-2,-1),3,2,byrow = T)
par(mfrow=c(1,3))
for (para in 1:3) {
  if (para==1) {v=c(0.01,0.05,0.1,0.5,1);xlab=expression(paste("True ",rho));ylab=expression(paste("Inferred ",rho))}
  if (para==2) {v=c(100,200,500,1000);xlab=expression(paste("True ",delta));ylab=expression(paste("Inferred ",delta))}
  if (para==3) {v=c(0.01,0.02,0.05,0.1);xlab=expression(paste("True ",nu));ylab=expression(paste("Inferred ",nu))}
  plot(c(),c(),xlim=ranges[para,],ylim=ranges[para,],xlab=xlab,ylab=ylab,xaxt="n",yaxt="n")
  v=log10(v)
  axis(1,at=v,labels=10^v)
  axis(2,at=v,labels=10^v)
  points(log10(allres[,(para-1)*7+1]),log10(allres[,(para-1)*7+5]),pch=20)
  lines(ranges[para,],ranges[para,])
  for (i in 1:nrow(allres))
    lines(c(log10(allres[i,(para-1)*7+1]),log10(allres[i,(para-1)*7+1])),c(log10(allres[i,(para-1)*7+6]),log10(allres[i,(para-1)*7+7])))
}
dev.off()
browseURL('ana_simu_multi.pdf')
