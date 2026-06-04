rm(list=ls())
library(RClonalFrame2)
#set.seed(1)
res=run_cf2(simu_pars=c(10,1e6,1e-4,0.1,1000,0.05,1),update_pars=c(0,0,0,0),burn=1000,iter=2000,thin=1,profile=T)
plot(res,'trace')
#plot(res,'simulation')

