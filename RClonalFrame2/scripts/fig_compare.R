library(RClonalFrame2)
rm(list=ls())
set.seed(0)
n=20
l=1000000
mu=0.001
niter=100
seed=round(runif(1,min=0,max=1e9))
system(sprintf('../build-release/comparison %d %d %f %d %d > res.txt',n,l,mu,niter,seed))
allres=read.table('res.txt',sep='\t',fill=T)
trees=allres[seq(1,nrow(allres),2),1]
allres=allres[seq(2,nrow(allres),2),]
system('rm res.txt')
t=read.tree(text=trees[1])
methods=c('pruning','joint','parsimony')
colnames(allres)<-methods
pdf('fig_compare.pdf',7,7)
plot(allres[,1],allres[,2],xlab='Unconditional likelihood',ylab='Conditional likelihood')
abline(0,1)
dev.off()
system('open fig_compare.pdf')
