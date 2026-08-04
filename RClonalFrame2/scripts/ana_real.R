rm(list=ls())
library(RClonalFrame2)
#alifile='~/data/st239/ali.fasta'
alifile='~/data/harris2010/core.full.aln'
load(sprintf('%s.RData',alifile))

plot(x[[3]],'trace',burnin=0)

library(coda)
m=as.mcmc(x[[3]])
print(effectiveSize(m))
gelman.diag(mcmc.list(as.mcmc(x[[3]]),as.mcmc(x[[4]])))

library(rwty)
r=as.rwty(x[[3]])
print(topological.approx.ess(r,burnin=0))

#Get dates from metadata
v=unlist(strsplit(alifile,'/'))
metafile=sprintf('%s/metadata.csv',paste0(v[1:(length(v)-1)],collapse = '/'))
tab=read.csv(metafile)
date=as.numeric(substr(tab$Date,1,4))#keep only year
names(date)=sprintf('%s-%s',tab$Genome.ID,tab$Genome.Name)

trees=list()
genome_length=length(x[[1]]$ancestral_sequences$patindex)
for (i in 1:4) {
  #Compute best tree for each run
  ind=medoid(x[[i]])
  tree=read.tree(text=x[[i]]$iterations$tree[ind])
  tree$edge.length=tree$edge.length*genome_length
  trees[[i]]=tree

  #Test temporal signal in each run
  library(BactDating)
  roottotip(initRoot(trees[[i]],date),date)
}
