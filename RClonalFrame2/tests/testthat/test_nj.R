#Testing basic package functions
context("Testing NJ algorithm")

test_that("Testing NJ algorithm.", {
  for (i in 1:10) {
    set.seed(i)
    tree=rtree(10)
    tree$edge.length=tree$edge.length*0.01
    ali=as.DNAbin(simSeq(tree,10000))
    write.FASTA(ali,'tmp.fasta')
    t2=nj(dist.dna(ali,model = 'JC69',pairwise.deletion = T))
    system('../../../build-debug/tests/test_nj')
    t3=read.tree('nj.nwk')
    expect_equal(path.dist(t2,t3,use.weight = T),0,tolerance=1e-4)
    system('rm tmp.fasta nj.nwk')
  }
})


