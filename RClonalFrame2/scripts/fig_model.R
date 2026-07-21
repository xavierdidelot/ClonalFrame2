library(RClonalFrame2)
rm(list=ls())
set.seed(1)
n=10#Number of sequences
l=10000#Genome length
rho=1
delta=1000
nu=0.01
b=0.001
x=run_cf2(simu_pars=c(n,l,rho,delta,nu,b,b,1),burn=0,iter=0,thin=1)
pdf('scripts/fig_model.pdf',10,5)
plot(x,'simulation_with_subs')
dev.off()
browseURL('scripts/fig_model.pdf')
