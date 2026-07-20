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
    if (tree$edge[i,2]<=n) bottom=seqs[[tree$tip.label[tree$edge[i,2]]]] else bottom=anc[[tree$node.label[tree$edge[i,2]-n]]]
    top=anc[[tree$node.label[tree$edge[i,1]-n]]]
    len=tree$edge.length[i]
    ex=exp(-len*4/3)
    pid =1/4+3/4*ex
    pnid=1/4-1/4*ex
    a=1-exp(-rho*len)
    #b=1-exp(-1/delta)
    if (a==0) b=0.5 else b=1-exp(-rho*len/expm1(rho*len*delta))#Using M/M/inf queue
    maxskip=100
    gamma=c(1-a,a,b,1-b)
    gammas=matrix(NA,maxskip,4)
    gammas[1,]=gamma
    if (maxskip>1) for (k in 2:maxskip) {
      gammas[k,1]=gammas[k-1,1]*pid*gamma[1]+gammas[k-1,2]*pid_nu*gamma[3]
      gammas[k,2]=gammas[k-1,1]*pid*gamma[2]+gammas[k-1,2]*pid_nu*gamma[4]
      gammas[k,3]=gammas[k-1,3]*pid*gamma[1]+gammas[k-1,4]*pid_nu*gamma[3]
      gammas[k,4]=gammas[k-1,3]*pid*gamma[2]+gammas[k-1,4]*pid_nu*gamma[4]
    }
    lscale=0

    #Forward algorithm
    is_sames=top[patindex]==bottom[patindex]
    jj=1
    todosites=sort(unique(c(seq(1,nsites,maxskip),nsites,which(!is_sames))))
    skippings=diff(c(0,todosites))
    lentodosites=length(todosites)
    alpha=matrix(NA,lentodosites+1,2)
    alpha[1,1]=b/(a+b)
    alpha[1,2]=a/(a+b)
    for (jj in 1:lentodosites) {
      j=todosites[jj]
      skipping=skippings[jj]
      if (is_sames[j]) {p1=pid;p2=pid_nu} else {p1=pnid;p2=pnid_nu}
      alpha[jj+1,1]=(alpha[jj,1]*gammas[skipping,1]+alpha[jj,2]*gammas[skipping,3])*p1
      alpha[jj+1,2]=(alpha[jj,1]*gammas[skipping,2]+alpha[jj,2]*gammas[skipping,4])*p2
      scale=alpha[jj+1,1]+alpha[jj+1,2]
      alpha[jj+1,]=alpha[jj+1,]/scale
      lscale=lscale+log(scale)
      j=j+1
    }
    ll=ll+lscale
  }
  ll=ll+nsites*log(1/4)
  return(ll)
}
