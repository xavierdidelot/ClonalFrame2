#Testing basic package functions
context("Testing pattern analysis")

test_that("Testing pattern analysis.", {
  for (i in 1:10) {
    set.seed(i)
    tree=rtree(10)
    tree$edge.length=tree$edge.length*0.01
    ali=as.DNAbin(simSeq(tree,10000))
    write.FASTA(ali,'tmp.fasta')
    m=system('../../../build-debug/tests/test_patterns',intern=T)
    m=as.numeric(m)
    pali=as.phyDat(ali)
    expect_equal(m,attr(pali,'nr'))
    system('rm tmp.fasta')
  }
})


