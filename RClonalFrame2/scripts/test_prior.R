# Testing that prior is sampled when no likelihood
rm(list=ls())
n=100#Number of sequences
topology_test=(n<10)#Slow if n>10
niter=1000#Number of iterations
bra_len_prior_mean=1#Mean of the exponential prior on branch lengths
x=run_cf2(simu_pars=c(n,0,bra_len_prior_mean,0,0,0,1),update_pars=c(100,0,0,0,0),burn=0,iter=niter,thin=1,profile=T,bra_len_prior_mean = bra_len_prior_mean)

trees=vector("list",niter)
b1=rep(0,niter)
for (i in 1:niter) {
  trees[[i]]=read.tree(text=x$iterations[i,'tree'])
  w=which(trees[[i]]$tip.label=='seq1')
  w=which(trees[[i]]$edge[,2]==w)
  b1[i]=trees[[i]]$edge.length[w]
}

if (topology_test) {
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
  if (length(counts)>1) print(chisq.test(counts))
}

trees2=vector("list",niter)
b2=rep(0,niter)
p2=rep(0,niter)
s2=rep(0,niter)
for (i in 1:niter) {
  trees2[[i]]=rtree(n,rooted=F,br=function(n) rexp(n,rate=1/bra_len_prior_mean))
  w=which(trees2[[i]]$tip.label=='t1')
  w=which(trees2[[i]]$edge[,2]==w)
  b2[i]=trees[[i]]$edge.length[w]
  s2[i]=sum(trees2[[i]]$edge.length)
  p2[i]=-s2[i]/bra_len_prior_mean-log(bra_len_prior_mean)*(2*n-3)
}

par(mfrow=c(3,3))
plot(x$iterations[,'tree_prior'],main='Tree prior',type='l',xlab='',ylab='')
plot(x$iterations[,'sum_bra_len'],main='Sum branch lengths',type='l',xlab='',ylab='')
plot(b1,main='Branch length 1',type='l',xlab='',ylab='')
hist(x$iterations[,'tree_prior'],main='Tree prior',xlab='',probability=T)
hist(x$iterations[,'sum_bra_len'],main='Sum branch lengths',xlab='',probability=T)
xs=seq(0,max(x$iterations[,'sum_bra_len']),length.out=1000)
lines(xs,dgamma(xs,shape=2*n-3,rate=1/bra_len_prior_mean))
hist(b1,main='Branch length 1',xlab='',probability=T)
xs=seq(0,max(b1),length.out=1000)
lines(xs,dexp(xs,rate=1/bra_len_prior_mean))

qqplot(p2,x$iterations[,'tree_prior'],xlab='rtree',ylab='CF2 prior')
abline(0, 1, col = 2, lty = 2)
qqplot(s2,x$iterations[,'sum_bra_len'],xlab='rtree',ylab='CF2 prior')
abline(0, 1, col = 2, lty = 2)
qqplot(b2,b1,xlab='rtree',ylab='CF2 prior')
abline(0, 1, col = 2, lty = 2)
