#Testing basic package functions
context("Testing parsimony algorithm")

test_that("Testing parsimony algorithm.", {
  for (i in 1:10) {
    set.seed(i)
    tree=rtree(10)
    tree$edge.length=tree$edge.length*0.01
    ali=as.DNAbin(simSeq(tree,10000))
    write.tree(tree,'tmp.nwk')
    write.FASTA(ali,'tmp.fasta')
    m=system('../../../build-debug/tests/test_parsimony',intern=T)
    m=as.numeric(unlist(strsplit(m,',')))
    expect_equal(m[1],m[2])
    pali=as.phyDat(ali)
    p=parsimony(tree, pali, method = "fitch")
    expect_equal(m[1],p)
    tree$node.label=paste0('N',1:tree$Nnode)
    anc=phangorn::ancestral.pars(tree,pali,return='ancestral')
    anc=anc[tree$node.label]
    r=likelihood_anc(tree,pali,anc,method='mine')
    print(c(m[3],r))
    #expect_equal(m[3],r)
    system('rm tmp.fasta tmp.nwk')
  }
})


