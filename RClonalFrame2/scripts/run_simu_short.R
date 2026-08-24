# Run on a simulated dataset
library(RClonalFrame2)
rm(list=ls())
set.seed(0)
n=20#Number of sequences
l=1e6#Genome length
rho=0.2
delta=1000
nu=0.05
x=run_cf2(simu_pars=c(n,l,rho,delta,nu,1e-4,1e-4),iter=1e3,profile=T,update_pars = c(1,1,1,1,1))
#plot(x)
#plot(x,'trace',burnin=0)
#plot(x,'trace_sd',burnin=0)
#plot(x,'trace_acceptance',burnin=0)
#plot(x,'trace_trees')
#plot(x,'simulation')
#print(tail(x$iterations$paccepted_spr,1))

filename=tempfile(fileext="pdf")
pdf(filename,14,7)
plot(x,'recombination')
dev.off()
browseURL(filename)
