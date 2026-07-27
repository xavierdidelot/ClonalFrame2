# Run on many simulated datasets
library(RClonalFrame2)
rm(list=ls())
n=100#Number of sequences
l=1e6#Genome length
for (ind in 1:100) {
  set.seed(ind)
  rho=runif(1,0.01,1)
  delta=runif(1,100,1000)
  nu=runif(1,0.01,0.1)
  set.seed(ind)
  x0=run_cf2(simu_pars=c(n,l,rho,delta,nu,1e-4,1e-4),iter=1e6,starting_pars=c(0,0,0,0),update_pars=c(0,0,0,1,1))
  set.seed(ind)
  x =run_cf2(simu_pars=c(n,l,rho,delta,nu,1e-4,1e-4),iter=1e6,starting_pars=c(0.1,500,0.05,0))
  save.image(sprintf('long%d.RData',ind))
}
