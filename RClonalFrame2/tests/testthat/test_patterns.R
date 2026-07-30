#Testing basic package functions
context("Testing pattern analysis")

test_that("Testing pattern analysis.", {
  for (i in 1:10) {
    l=100
    set.seed(i)
    tree=rtree(10)
    tree$edge.length=tree$edge.length*0.0#0.01
    ali=as.DNAbin(simSeq(tree,l))
    for (i in seq_len(nrow(ali))) {
      idx <- sample(ncol(ali), round(ncol(ali) * 0.5))
      ali[i, idx] <- as.raw(0xf0)
    }
    write.FASTA(ali,'tmp.fasta')
    #system(sprintf('echo ">empty\n%s" >> tmp.fasta',paste0(rep('N',l),collapse = '')))
    m=system('../../../build-debug/test_patterns',intern=T)
    m=as.numeric(m)
    ali=read.FASTA('tmp.fasta')
    pali=as.phyDat(ali)
    expect_equal(m,attr(pali,'nr'))
    #system('rm tmp.fasta')
  }
})


