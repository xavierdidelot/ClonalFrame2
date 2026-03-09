#Testing basic package functions
context("Testing jointml algorithm")

test_that("Testing jointml algorithm.", {
  for (i in 1:10) {
    set.seed(i)
    tree=rtree(10)
    tree$edge.length=tree$edge.length*0.01
    ali=as.DNAbin(simSeq(tree,10000))
    write.tree(tree,'tmp.nwk')
    write.FASTA(ali,'tmp.fasta')
    m=system('../../../build-debug/tests/test_jointml',intern=T)
    m=as.numeric(unlist(strsplit(m,',')))
    expect_equal(m[1],m[2])
    pali=as.phyDat(ali)
    p=likelihood(tree,pali,method='joint')
    expect_equal(m[1],p)
    system('rm tmp.fasta tmp.nwk')
  }
})


