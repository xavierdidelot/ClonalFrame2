# Analysis of many simulated datasets
library(RClonalFrame2)
rm(list=ls())
n=100#Number of sequences
l=1e6#Genome length
burnin=0.5
allres=c()
nams=c('rho','delta','nu')
scores=c()
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

  truetree=read.tree(text=x0$simulation_parameters$tree)
  l0=length(x0$iterations$tree)
  l1=length(x$iterations$tree)
  tree0=read.tree(text=x0$iterations$tree[l0],1)
  tree1=read.tree(text= x$iterations$tree[l1],1)
  score0=KF.dist(truetree,tree0)
  score1=KF.dist(truetree,tree1)
  score0int=KF.dist(tree0,read.tree(text=x0$iterations$tree[l0/2],1))
  score1int=KF.dist(tree1,read.tree(text= x$iterations$tree[l1/2],1))
  tree0$edge.length=tree0$edge.length/sum(tree0$edge.length)*sum(truetree$edge.length)
  tree1$edge.length=tree1$edge.length/sum(tree1$edge.length)*sum(truetree$edge.length)
  score0scaled=KF.dist(truetree,tree0)
  score1scaled=KF.dist(truetree,tree1)
  scores=rbind(scores,c(score0,score1,score0int,score1int,score0scaled,score1scaled))
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

para=1
length(which(allres[,(para-1)*7+1]>=allres[,(para-1)*7+6] & allres[,(para-1)*7+1]<=allres[,(para-1)*7+7]))

colMeans(scores)
