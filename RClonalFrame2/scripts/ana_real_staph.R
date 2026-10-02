rm(list=ls())
library(BactDating)
library(RClonalFrame2)
#alifile='~/data/st239/ali.fasta'
alifile='~/data/harris2010/core.full.aln'
load(sprintf('%s.RData',alifile))

#Make figure
pdf(sprintf('%s-fig.pdf',alifile),10,10)
par(fig=c(0,0.5,0,0.5),mai=c(0.5,0.5,0.5,0.5))
roottotip(initRoot(trees[[4]],date),date,showText = F,showTree = F,xl=c(1960,2007))
par(fig=c(0.45,1,0,0.5),new=T)
plot(bd[[4]],'treeCI',show.tip.label=F,xl=c(1940,2007))
par(fig=c(0,0.5,0.5,1),new=T)
roottotip(initRoot(trees[[1]],date),date,showText = F,showTree = F,xl=c(1960,2007))
par(fig=c(0.45,1,0.5,1),new=T)
plot(bd[[1]],'treeCI',show.tip.label=F,xl=c(1940,2007))
par(xpd=NA)
text(-0.1,1.05,'A',cex=2)
text(0.5,1.05,'B',cex=2)
text(-0.1,0.45,'C',cex=2)
text(0.5,0.45,'D',cex=2)
dev.off()
system(sprintf('open %s-fig.pdf',alifile))

