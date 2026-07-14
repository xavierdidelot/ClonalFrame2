# Run on a simulated dataset
library(RClonalFrame2)
rm(list=ls())
set.seed(0)
n=100#Number of sequences
l=1e6#Genome length
rho=0.1
delta=500
nu=0.05
x=run_cf2(simu_pars=c(n,l,rho,delta,nu,1e-4,1e-4),iter=1e5,profile=T)
save.image('long.RData')
plot(x)
plot(x,'trace',burnin=0)
plot(x,'trace_sd',burnin=0)
plot(x,'trace_acceptance',burnin=0)
#plot(x,'trace_trees')
#plot(x,'simulation')

library(coda)
m=as.mcmc(x)
print(effectiveSize(m))

library(rwty)
r=as.rwty(x)
print(topological.approx.ess(r,burnin=0))

m=as.multiPhylo(x)
ind=medoid(m)
cons_topo=consensus(m, p = 1, rooted = FALSE)
cons=add_edge_length(cons_topo, m, rooted = FALSE)
plot(cons,show.tip.label = F)
axisPhylo(1)
