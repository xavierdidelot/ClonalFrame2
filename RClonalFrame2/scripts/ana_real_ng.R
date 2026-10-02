library(ape)
if (F) {
rm(list=ls())
library(RClonalFrame2)
alifile='~/data/grad2014/clade1.fasta'
load(sprintf('%s.RData',alifile))
allx=x
x=x[[3]]
iteration=medoid(x)
r=decode_recombination(x,iteration)
}

#pdf(sprintf('%s.pdf',alifile),20,10)
png(sprintf('%s.png',alifile), width = 5000, height = 3000, res = 300)
newick=x$iterations$tree[iteration]
#newick=write.tree(unroot(phangorn::midpoint(read.tree(text=newick))))
genome_length=length(x$ancestral_sequences$patindex)
tree_aligned=RClonalFrame2:::plot_tree_aligned(newick,genome_length,show.tip.label = F)

#Reorder rows of matrix r to match order of leaves in tree_aligned
treeit=read.tree(text=newick)
treeit_strings=RClonalFrame2:::children_strings(treeit)
r2=matrix(0,nrow(r),ncol(r))
for (i in 1:Ntip(tree_aligned)) {
  str=tree_aligned$tip.label_strings[i]
  w=which(treeit_strings==str)
  if (length(w)==0) next
  w=which(treeit$edge[,2]==w)
  r2[i,]=r[w,]
}

#r2=r2[,1:20000]
subint=5
r3=matrix(NA,nrow(r2),ceiling(ncol(r2)/subint))
for (i in 1:ncol(r3)) {
  sub=r2[,max(1,(i-3)*subint):min((i+3)*subint,ncol(r2))]
  r3[,i]=apply(sub,1,max)
}
image(x=seq(1,genome_length,length.out=ncol(r3)),y=1:nrow(r3),z=t(r3),useRaster=F,axes = F,add=T)

dev.off()
system(sprintf('open %s.png',alifile))

