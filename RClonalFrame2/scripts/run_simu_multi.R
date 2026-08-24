# Run on many simulated datasets
library(RClonalFrame2)
rm(list=ls())
n=100#Number of sequences
l=1e6#Genome length
for (ind in 1:100) {
  set.seed(ind)
  rho=10^runif(1,-2,0)
  delta=10^runif(1,2,3)
  nu=10^runif(1,-2,-1)
  set.seed(ind)
  x0=run_cf2(simu_pars=c(n,l,rho,delta,nu,1e-4,1e-4),iter=1e5,starting_pars=c(0,0,0),update_pars=c(0,0,0,1,1))
  set.seed(ind)
  x =run_cf2(simu_pars=c(n,l,rho,delta,nu,1e-4,1e-4),iter=1e5,starting_pars=c(0.1,500,0.05))
  save.image(sprintf('long%d.RData',ind))
}
