#Testing basic package functions
context("Testing pattern analysis")

test_that("Testing pattern analysis.", {
  source('utils_testing.R')
  for (rep in 1:10) {
    set.seed(rep)
    tree=rtree(10)
    tree$edge.length=tree$edge.length*0.1
    ali=as.DNAbin(simSeq(tree,1000))
    for (i in seq_len(nrow(ali))) {
      idx <- sample(ncol(ali), round(ncol(ali) * 0.5))
      ali[i, idx] <- as.raw(0xf0)
    }
    tdir=tempdir()
    write.FASTA(ali,sprintf('%s/tmp.fasta',tdir))
    exe=normalizePath('../../../build-debug/test_patterns')
    m=system(sprintf('cd %s;%s',tdir,exe),intern=T)
    m=as.numeric(m)
    ali=read.FASTA(sprintf('%s/tmp.fasta',tdir))
    pali=dnabin_to_phydat_collapsing(ali)
    expect_equal(m,attr(pali,'nr'))
  }
})


