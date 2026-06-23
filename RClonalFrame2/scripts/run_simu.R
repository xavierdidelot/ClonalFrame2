# Run on a simulated dataset
library(RClonalFrame2)
rm(list=ls())
set.seed(0)
n=50#Number of sequences
l=1e6#Genome length
rho=0
delta=500
nu=0.05
bra_len_prior_mean=1e-4#Mean of the exponential prior on branch lengths
x=run_cf2(simu_pars=c(n,l,bra_len_prior_mean,rho,delta,nu,1),update_pars=c(10,0,0,0,10),burn=100,iter=200,thin=1,profile=T,bra_len_prior_mean = bra_len_prior_mean)
plot(x)
plot(x,'trace_acceptance')
#plot(x,'trace_trees')
#plot(x,'simulation')

#final acceptance value:
tail(x$iterations[,'paccepted_tree'],1)
