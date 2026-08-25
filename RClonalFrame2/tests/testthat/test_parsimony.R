#Testing basic package functions
context("Testing parsimony algorithm")

test_that("Testing parsimony algorithm.", {
  for (i in 1:10) {
    set.seed(i)
    tree=rtree(10)
    tree$edge.length=tree$edge.length*0.01
    ali=as.DNAbin(simSeq(tree,10000))
    tdir=tempdir()
    write.tree(tree,sprintf('%s/tmp.nwk',tdir))
    write.FASTA(ali,sprintf('%s/tmp.fasta',tdir))
    exe=normalizePath('../../../build-debug/test_parsimony')
    m=system(sprintf('cd %s;%s',tdir,exe),intern=T)
    m=as.numeric(unlist(strsplit(m,',')))
    expect_equal(m[1],m[2])
    pali=as.phyDat(ali)
    p=parsimony(tree, pali, method = "fitch")
    expect_equal(m[1],p)
  }
})


