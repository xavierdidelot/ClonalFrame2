#Testing basic package functions
context("Testing likelihood functions")

test_that("Testing likelihood functions.", {
  for (i in 1:10) {
    set.seed(i)
    tree=rtree(10)
    tree$edge.length=tree$edge.length*0.01
    ali=as.DNAbin(simSeq(tree,10000))
    write.tree(tree,'tmp.nwk')
    write.FASTA(ali,'tmp.fasta')
    m=system('../../../build-debug/tests/test_likelihoods',intern=T)
    a=as.numeric(unlist(strsplit(m,',')))
    tree2=read.tree('tmp2.nwk')
    anc2=as.phyDat(read.FASTA('tmp2.fasta'))
    seqs2=anc2[tree2$tip.label]
    anc2=anc2[tree2$node.label]
    b=likelihood_anc(tree2,seqs2,anc2,method='mine')
    expect_equal(a,b)
    system('rm tmp.fasta tmp.nwk tmp2.fasta tmp2.nwk')
  }
})
