#Testing basic package functions
context("Testing likelihood functions")

test_that("Testing likelihood functions.", {
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
    if (i%%3==2) tree=unroot(tree)
    if (i%%3==0) tree=midpoint(tree)
    write.tree(tree,'tmp.nwk')
    write.FASTA(ali,'tmp.fasta')
    m=system('../../../build-debug/test_likelihoods',intern=T)
    a=as.numeric(unlist(strsplit(m,',')))
    pali=dnabin_to_phydat_collapsing(ali)
    ali=as.DNAbin(pali)
    ll=likelihood(tree,pali,method='pruning')
    ll2=likelihood(midpoint(tree),pali,method='pruning')
    expect_equal(ll,ll2)
    expect_equal(a[1],ll)

    tree2=read.tree('tmp2.nwk')
    expect_equal(path.dist(tree,tree2,use.weight = T),0,tolerance=1e-6)
    anc2=as.phyDat(read.FASTA('tmp2.fasta'))#contains uncertainties
    tmp=as.DNAbin(anc2)
    for (i in tree$tip.label) expect_equal(all(tmp[i,]==ali[i,]),T)
    anc2=as.phyDat(read.FASTA('tmp3.fasta'))#does not contain uncertainties
    seqs2=anc2[tree2$tip.label]
    anc2=anc2[tree2$node.label]
    b=likelihood_anc(tree2,seqs2,anc2,method='mine')
    c=likelihood_anc(tree2,seqs2,anc2,method='phangorn')
    expect_equal(b,c)
    expect_equal(a[2],b)
    d=likelihood_anc_rec(tree2,seqs2,anc2,rho=0,delta=202,nu=0.055)
    expect_equal(b,d)
    expect_equal(a[3],b)
    e=likelihood_anc_rec(tree2,seqs2,anc2,rho=1.1,delta=202,nu=0.055)
    expect_equal(a[4],e)
    system('rm tmp.fasta tmp.nwk tmp2.fasta tmp2.nwk tmp3.fasta')
  }
})
