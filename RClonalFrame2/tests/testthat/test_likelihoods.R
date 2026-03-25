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
    ll=likelihood(tree,as.phyDat(ali),method='pruning')
    ll2=likelihood(midpoint(tree),as.phyDat(ali),method='pruning')
    expect_equal(ll,ll2)
    expect_equal(a[1],ll)

    tree2=read.tree('tmp2.nwk')
    expect_equal(path.dist(tree,tree2,use.weight = T),0,tolerance=1e-6)
    anc2=as.phyDat(read.FASTA('tmp2.fasta'))
    tmp=as.DNAbin(anc2)
    for (i in tree$tip.label) expect_equal(all(tmp[i,]==ali[i,]),T)
    seqs2=anc2[tree2$tip.label]
    anc2=anc2[tree2$node.label]
    b=likelihood_anc(tree2,seqs2,anc2,method='mine')
    c=likelihood_anc(tree2,seqs2,anc2,method='phangorn')
    expect_equal(b,c)
    expect_equal(a[2],b)
    system('rm tmp.fasta tmp.nwk tmp2.fasta tmp2.nwk')
  }
})
