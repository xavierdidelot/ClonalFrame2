library(RClonalFrame2)
#alifile='~/data/st239/ali.fasta'
alifile='~/data/harris2010/core.full.aln'
iter=1e5
x0=run_cf2(alifile,iter=iter,starting_pars=c(0,0,0,0),update_pars=c(0,0,0,1,1),seed=0)
x1=run_cf2(alifile,iter=iter,starting_pars=c(0,0,0,0),update_pars=c(0,0,0,1,1),seed=1)
x2=run_cf2(alifile,iter=iter,seed=2)
x2=run_cf2(alifile,iter=iter,seed=3)
save.image(sprintf('%s.Rdata',alifile))

#plot(x)
#plot(x,'trace')
#plot(x,'trace_acceptance')
