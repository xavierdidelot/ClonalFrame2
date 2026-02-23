#Testing basic package functions
context("Testing likelihood functions")

test_that("Testing likelihood functions.", {
  set.seed(0)
  tree=rtree(10)
  tree$edge.length=tree$edge.length*0.01
  seqs=phangorn::simSeq(tree,l=100,rate=0)
  tree$node.label=paste0('N',1:tree$Nnode)
  anc=phangorn::ancestral.pars(tree,seqs,return='ancestral')
  anc=anc[tree$node.label]
  #anc=phangorn:::joint_pml(phangorn::pml(tree,seqs,model='JC'))
  a=likelihood_anc(tree,seqs,anc,method='phangorn')
  b=likelihood_anc(tree,seqs,anc,method='mine')
  c=100*(sum(log(1/4+3/4*exp(-tree$edge.length)))+log(1/4))
  expect_equal(b,c)#but a is not the same...
})
