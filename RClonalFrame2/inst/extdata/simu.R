#Script used to create the example data files
rm(list=ls())
set.seed(0)
library(ape)
library(phangorn)
tree=rtree(20)
mu=0.01#mutation rate per site per branch length
tree$edge.length=tree$edge.length*mu
write.tree(tree,'inst/extdata/simu.nwk')
seqs=simSeq(tree,l=10000)
write.phyDat(seqs,'inst/extdata/simu.fas',format="fasta")
save(tree,seqs,file='data/simu.rda')
