# Testing that prior is sampled when no likelihood
rm(list=ls())
library(RClonalFrame2)
set.seed(1)
n=100#Number of sequences
topology_test=(n<10)#Slow if n>10
niter=2000#Number of iterations
prior_pars=c(1,0.5,500,250,0.05,0.025,0.001,0.0005)
x=run_cf2(simu_pars=c(n,0,prior_pars[1],prior_pars[3],prior_pars[5],prior_pars[7],prior_pars[8]),update_pars=c(100,100,100,100,100),burn=0,iter=niter,thin=1,prior_pars = prior_pars)

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

b_shape=prior_pars[7]^2/prior_pars[8]^2
b_scale=prior_pars[8]^2/prior_pars[7]
trees2=vector("list",niter)
b2=rep(0,niter)
p2=rep(0,niter)
s2=rep(0,niter)
for (i in 1:niter) {
  trees2[[i]]=rtree(n,rooted=F,br=function(n) rgamma(n,shape=b_shape,scale=b_scale))
  w=which(trees2[[i]]$tip.label=='t1')
  w=which(trees2[[i]]$edge[,2]==w)
  b2[i]=trees[[i]]$edge.length[w]
  s2[i]=sum(trees2[[i]]$edge.length)
  p2[i]=sum(dgamma(trees2[[i]]$edge.length,shape=b_shape,scale=b_scale,log=T))
}

source('scripts/fig_label.R')
pdf('scripts/test_prior.pdf',10,15)
par(mfrow=c(6,3))
plot(x$iterations[,'prior_tree'],main='Tree prior',type='l',xlab='',ylab='')
fig_label('A',cex=2)
hist(x$iterations[,'prior_tree'],main='Tree prior',xlab='',probability=T,breaks=20)
fig_label('B',cex=2)
qqplot(p2,x$iterations[,'prior_tree'],xlab='Expectation',ylab='Samples',main='Tree prior')
fig_label('C',cex=2)
abline(0, 1, col = 2, lty = 2)

plot(x$iterations[,'sum_bra_len'],main='Sum branch lengths',type='l',xlab='',ylab='')
fig_label('D',cex=2)
hist(x$iterations[,'sum_bra_len'],main='Sum branch lengths',xlab='',probability=T,breaks=20)
fig_label('E',cex=2)
xs=seq(min(x$iterations[,'sum_bra_len']),max(x$iterations[,'sum_bra_len']),length.out=1000)
lines(xs,dgamma(xs,shape=(2*n-3)*b_shape,scale=b_scale))
qqplot(s2,x$iterations[,'sum_bra_len'],xlab='Expectation',ylab='Samples',main='Sum branch lengths')
fig_label('F',cex=2)
abline(0, 1, col = 2, lty = 2)

plot(b1,main='Branch length 1',type='l',xlab='',ylab='')
fig_label('G',cex=2)
hist(b1,main='Branch length 1',xlab='',probability=T,breaks=20)
fig_label('H',cex=2)
xs=seq(min(b1),max(b1),length.out=1000)
lines(xs,dgamma(xs,shape=b_shape,scale=b_scale))
qqplot(b2,b1,xlab='Expectation',ylab='Samples',main='Branch length 1')
fig_label('I',cex=2)
abline(0, 1, col = 2, lty = 2)

rho_shape=prior_pars[1]^2/prior_pars[2]^2
rho_scale=prior_pars[2]^2/prior_pars[1]
delta_shape=prior_pars[3]^2/prior_pars[4]^2
delta_scale=prior_pars[4]^2/prior_pars[3]
nu_shape=prior_pars[5]^2/prior_pars[6]^2
nu_scale=prior_pars[6]^2/prior_pars[5]

plot(x$iterations[,'rho'],main='Rho',type='l',xlab='',ylab='')
fig_label('J',cex=2)
hist(x$iterations[,'rho'],main='Rho',xlab='',probability=T,breaks=20)
fig_label('K',cex=2)
xs=seq(0,max(x$iterations[,'rho']),length.out=1000)
lines(xs,dgamma(xs,shape=rho_shape,scale=rho_scale))
qqplot(rgamma(nrow(x$iterations),shape=rho_shape,scale=rho_scale),x$iterations[,'rho'],xlab='Expectation',ylab='Samples',main='Rho')
fig_label('L',cex=2)
abline(0, 1, col = 2, lty = 2)

plot(x$iterations[,'delta'],main='Delta',type='l',xlab='',ylab='')
fig_label('M',cex=2)
hist(x$iterations[,'delta'],main='Delta',xlab='',probability=T,breaks=20)
fig_label('N',cex=2)
xs=seq(0,max(x$iterations[,'delta']),length.out=1000)
lines(xs,dgamma(xs,shape=delta_shape,scale=delta_scale))
qqplot(rgamma(nrow(x$iterations),shape=delta_shape,scale=delta_scale),x$iterations[,'delta'],xlab='Expectation',ylab='Samples',main='Delta')
abline(0, 1, col = 2, lty = 2)
fig_label('O',cex=2)

plot(x$iterations[,'nu'],main='Nu',type='l',xlab='',ylab='')
fig_label('P',cex=2)
hist(x$iterations[,'nu'],main='Nu',xlab='',probability=T,breaks=20)
fig_label('Q',cex=2)
xs=seq(0,max(x$iterations[,'nu']),length.out=1000)
lines(xs,dgamma(xs,shape=nu_shape,scale=nu_scale))
qqplot(rgamma(nrow(x$iterations),shape=nu_shape,scale=nu_scale),x$iterations[,'nu'],xlab='Expectation',ylab='Samples',main='Nu')
abline(0, 1, col = 2, lty = 2)
fig_label('R',cex=2)
dev.off()
system('open scripts/test_prior.pdf')
