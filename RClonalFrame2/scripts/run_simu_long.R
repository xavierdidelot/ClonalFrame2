# Run on a simulated dataset
library(RClonalFrame2)
rm(list=ls())
set.seed(0)
n=50#Number of sequences
l=1e6#Genome length
rho=0
delta=500
nu=0.05
x=run_cf2(simu_pars=c(n,l,rho,delta,nu,1e-4,1e-4),burn=1e5,iter=2e5,thin=1e2,profile=T)
plot(x)
plot(x,'trace_acceptance')
#plot(x,'trace_trees')
#plot(x,'simulation')
