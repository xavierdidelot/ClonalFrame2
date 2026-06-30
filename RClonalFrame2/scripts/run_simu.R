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
x=run_cf2(simu_pars=c(n,l,rho,delta,nu,bra_len_prior_mean,bra_len_prior_mean),update_pars=c(0,0,0,10,10),burn=1000,iter=2000,thin=1,profile=T)
plot(x)
plot(x,'trace_acceptance')
#plot(x,'trace_trees')
#plot(x,'simulation')
