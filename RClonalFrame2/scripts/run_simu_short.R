# Run on a simulated dataset
library(RClonalFrame2)
rm(list=ls())
set.seed(0)
n=5#Number of sequences
l=1e6#Genome length
rho=1
delta=500
nu=0.05
x=run_cf2(simu_pars=c(n,l,rho,delta,nu,1e-4,1e-4),iter=1e3,profile=T)
plot(x)
plot(x,'trace_acceptance')
#plot(x,'trace_trees')
#plot(x,'simulation')
