#' Decoding of recombined regions using forward-backward algorithm
#' @param tree Phylogeny in ape::phylo format
#' @param seqs Sequences in phangorn::phyDat format
#' @param anc Ancestral sequences in phangorn::phyDat format
#' @param rho Rho parameter
#' @param delta Delta parameter
#' @param nu Nu parameter
#' @returns Probability of recombination for each branch and site
#'
decode_recombination <- function(tree,seqs,anc,rho,delta,nu) {
  ll=0
  n=Ntip(tree)
  w=attr(seqs,'weight')
  patindex=attr(seqs,'index')
  nsites=length(patindex)
  pid_nu =1/4+3/4*exp(-nu*4/3)
  pnid_nu=1/4-1/4*exp(-nu*4/3)
  for (i in 1:nrow(tree$edge)) {
    if (tree$edge[i,2]<=n) bottom=seqs[[tree$tip.label[tree$edge[i,2]]]]
    else bottom=anc[[tree$node.label[tree$edge[i,2]-n]]]
    top=anc[[tree$node.label[tree$edge[i,1]-n]]]
    len=tree$edge.length[i]
    ex=exp(-len*4/3)
    pid =1/4+3/4*ex
    pnid=1/4-1/4*ex
    a=1-exp(-rho*len)
    #b=1-exp(-1/delta)
    if (a==0) b=0.5 else b=1-exp(-rho*len/expm1(rho*len*delta))#Using M/M/inf queue
    alpha1=b/(a+b)
    alpha2=a/(a+b)
    lscale=0
    #Forward algorithm
    for (j in 1:nsites) {
      pat=patindex[j]
      is_same=top[pat]==bottom[pat]
      if (is_same) {p1=pid;p2=pid_nu} else {p1=pnid;p2=pnid_nu}
      newalpha1=(alpha1*(1-a)+alpha2*b)*p1
      alpha2=(alpha1*a+alpha2*(1-b))*p2
      alpha1=newalpha1
      scale=alpha1+alpha2
      alpha1=alpha1/scale
      alpha2=alpha2/scale
      lscale=lscale+log(scale)
    }
    ll=ll+lscale
  }
  ll=ll+nsites*log(1/4)
  return(ll)
}
