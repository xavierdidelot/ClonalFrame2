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
if (alifile=='~/data/harris2010/core.full.aln') date=c(date,Reference=2003)#TW20 reference

trees=list()
bd=list()
genome_length=length(x[[1]]$ancestral_sequences$patindex)
for (i in 1:4) {
  #Compute best tree for each run
  ind=medoid(x[[i]])
  tree=read.tree(text=x[[i]]$iterations$tree[ind])
  tree$edge.length=tree$edge.length*genome_length
  trees[[i]]=tree
}

#pdf(sprintf('%s.pdf',alifile),10,8)
#plot(x[[3]],'recombination')
#dev.off()

library(BactDating)
for (i in 1:4) {
  #Test temporal signal in each run
  roottotip(initRoot(trees[[i]],date),date)
}

for (i in 1:4) {
  #Run bactdate for each run
  bd[[i]]=bactdate(trees[[i]],date,showProgress=T,nbIts=1e6)
}

library(DiagnoDating)
for (i in 1:4) {
  print(ppcheck(bd[[i]],showPlot = F,showProgress = F))
  postdistpvals(bd[[i]],showPlot = F)
}

#Make figure
pdf(sprintf('%s-fig.pdf',alifile),10,10)
par(fig=c(0,0.5,0,0.5),mai=c(0.5,0.5,0.5,0.5))
roottotip(initRoot(trees[[4]],date),date,showText = F,showTree = F)
par(fig=c(0.5,1,0,0.5),new=T)
plot(bd[[4]],'treeCI',show.tip.label=F)
par(fig=c(0,0.5,0.5,1),new=T)
roottotip(initRoot(trees[[1]],date),date,showText = F,showTree = F)
par(fig=c(0.5,1,0.5,1),new=T)
plot(bd[[1]],'treeCI',show.tip.label=F)
par(xpd=NA)
text(-0.1,1.05,'A',cex=2)
text(0.5,1.05,'B',cex=2)
text(-0.1,0.45,'C',cex=2)
text(0.5,0.45,'D',cex=2)
dev.off()
system(sprintf('open %s-fig.pdf',alifile))

