#' Decoding of recombined regions using forward-backward algorithm
#' @param x Output from ClonalFrame2
#' @param iteration Index of the iteration to extract (default is the last iteration)
#' @param gap Gap between sites in output (default is 100)
#' @returns Matrix of probability of recombination for each branch and n-th site with n=gap
#' @export
#'
decode_recombination <- function(x,iteration,gap) {
  stopifnot(inherits(x, "cf2out"))
  if (missing(iteration)) iteration=nrow(x$iterations)
  if (missing(gap)) gap=100

  #Extract properties of iteration
  rho=x$iterations[iteration,'rho']
  delta=x$iterations[iteration,'delta']
  nu=x$iterations[iteration,'nu']
  tree=read.tree(text=x$iterations[iteration,'tree'])
  likelihood=x$iterations[iteration,'likelihood']

  #Recover data
  ancseqs=x$ancestral_sequences
  a=ancseqs$patmatrix$data[1:Ntip(tree)]
  names(a)=ancseqs$patmatrix$branch[1:Ntip(tree)]
  l=as.list(a)
  seqs_split=lapply(l, function(x) strsplit(x, "")[[1]])
  seqs=phyDat(seqs_split, type = "DNA")
  tmpindex=attr(seqs,'index')
  attr(seqs,'index')=tmpindex[ancseqs$patindex+1]
  attr(seqs,'weight')=tabulate(attr(seqs,'index'), nbins = attr(seqs,'nr'))

  #Perform ancestral state reconstruction
  pml=phangorn::pml(tree,seqs,model='JC')
  anc=phangorn:::joint_pml(pml)

  #Impute missing data if needed
  #seqs=resolve_tip_uncertainties(tree,seqs,anc)

  #Inference of recombined regions
  ll=0
  ll2=0
  n=Ntip(tree)
  w=attr(seqs,'weight')
  patindex=attr(seqs,'index')
  nsites=length(patindex)
  pid_nu =1/4+3/4*exp(-nu*4/3)
  pnid_nu=1/4-1/4*exp(-nu*4/3)
  outsites=seq(1,nsites,gap)
  res=matrix(NA,nrow(tree$edge),length(outsites))
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

    #Forward algorithm
    lscale=0
    is_sames=top[patindex]==bottom[patindex]
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
      scale=sum(alpha[jj+1,])
      alpha[jj+1,]=alpha[jj+1,]/scale
      lscale=lscale+log(scale)
    }
    ll=ll+lscale

    #Backward algorithm
    lscale2=0
    beta=c(1,1)
    for (jj in (lentodosites-1):1) {
      j=todosites[jj]
      skipping=skippings[jj+1]
      if (is_sames[todosites[jj+1]]) {p1=pid;p2=pid_nu} else {p1=pnid;p2=pnid_nu}
      newbeta1=beta[1]*gammas[skipping,1]*p1+beta[2]*gammas[skipping,2]*p2
      beta[2] =beta[1]*gammas[skipping,3]*p1+beta[2]*gammas[skipping,4]*p2
      beta[1]=newbeta1
      scale=sum(beta)
      beta=beta/scale
      lscale2=lscale2+log(scale)
      alpha[jj+1,1]=alpha[jj+1,1]*beta[1]
      alpha[jj+1,2]=alpha[jj+1,2]*beta[2]
      alpha[jj+1,]=alpha[jj+1,]/sum(alpha[jj+1,])
    }
    if (is_sames[todosites[1]]) {p1=pid;p2=pid_nu} else {p1=pnid;p2=pnid_nu}
    ll2 <- ll2 + log(b/(a+b)*beta[1]*p1 + a/(a+b)*beta[2]*p2) + lscale2

    #Interpolation for output
    res[i,]=approx(todosites,alpha[,2][-1],outsites)$y
  }
  ll=ll+nsites*log(1/4)
  ll2=ll2+nsites*log(1/4)

  #Check that all four likelihoods agree
  #d=likelihood_anc_rec(tree,seqs,anc,rho=rho,delta=delta,nu=nu)#slow
  #print(c(likelihood,d,ll,ll2))

  return(res)
}
