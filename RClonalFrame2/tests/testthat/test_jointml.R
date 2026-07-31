#Testing basic package functions
context("Testing jointml algorithm")

test_that("Testing jointml algorithm.", {
  source('utils_testing.R')
  for (i in 1:10) {
    set.seed(i)
    tree=rtree(10)
    tree$edge.length=tree$edge.length*0.01
    ali=as.DNAbin(simSeq(tree,10000))
    for (i in seq_len(nrow(ali))) {
      idx <- sample(ncol(ali), round(ncol(ali) * 0.5))
      ali[i, idx] <- as.raw(0xf0)
    }
    write.tree(tree,'tmp.nwk')
    write.FASTA(ali,'tmp.fasta')
    m=system('../../../build-debug/test_jointml',intern=T)
    m=as.numeric(unlist(strsplit(m,',')))
    expect_equal(m[1],m[2])
    pali=as.phyDat(ali)
    #pali=dnabin_to_phydat_collapsing(ali)
    p=likelihood(tree,pali,method='joint')
    p2=likelihood(midpoint(tree),pali,method='joint')
    expect_equal(p,p2)
    expect_equal(m[1],p)
    system('rm tmp.fasta tmp.nwk')
  }
})


