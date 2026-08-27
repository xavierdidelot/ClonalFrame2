#' Load sequences from FASTA file
#' @param filename Name of input FASTA file
#' @returns Sequences in phyDat format
#'
loadSeqs <- function(filename) {
  x=phangorn::read.phyDat(filename,format='fasta')
  return(x)
}

#' Likelihood function given ancestral sequences
#' @param tree Phylogeny in ape::phylo format
#' @param seqs Sequences in phangorn::phyDat format
#' @param anc Ancestral sequences in phangorn::phyDat format
#' @param workenv Work environment, optional but useful to speed up
#' @param method Method to use for calculation, can be either phangorn or mine
#' @returns Log-likelihood
#'
likelihood_anc <- function(tree,seqs,anc,workenv=NULL,method='mine') {
  seqs=resolve_tip_uncertainties(tree,seqs,anc)
  if (method=='phangorn') {
    nanc=names(anc)#setdiff(names(anc),tree$node.label[1])#remove root
    if (is.null(workenv)) workenv=list()#use a list, not an environment
    new_tree=phangorn::add.tips(tree=tree,tips=nanc,where=nanc,edge.length=rep(0,length(anc)))
    #both=phangorn:::rbind.phyDat(seqs,anc)#slow so replaced with below which assumes same patterns
    both=seqs;both[nanc]=anc[nanc]
    if (is.null(workenv$pml2)) {
      workenv$pml2=phangorn::pml(new_tree,both,model='JC')
    } else {
      workenv$pml2=stats::update(workenv$pml2,tree=new_tree,data=both)
    }
    return(workenv$pml2$logLik)
  }

  if (method=='mine') {
    ll=0
    n=Ntip(tree)
    w=attr(seqs,'weight')
    nsites=length(attr(seqs,'index'))
    for (i in 1:nrow(tree$edge)) {
      if (tree$edge[i,2]<=n) bottom=seqs[[tree$tip.label[tree$edge[i,2]]]]
      else bottom=anc[[tree$node.label[tree$edge[i,2]-n]]]
      top=anc[[tree$node.label[tree$edge[i,1]-n]]]
      l=tree$edge.length[i]
      ex=exp(-l*4/3)
      pid =log(1/4+3/4*ex)
      pnid=log(1/4-1/4*ex)
      ll=ll+pid *sum(w[top==bottom])
      ll=ll+pnid*sum(w[top!=bottom])
    }
    ll=ll+nsites*log(1/4)
    return(ll)
  }
}

#' Likelihood function
#' @param tree Phylogeny in ape::phylo format
#' @param seqs Sequences in phangorn::phyDat format
#' @param method Which method to use: pruning, joint, marginal or parsimony
#' @param workenv Work environment, optional but useful to speed up
#' @returns Log-likelihood
#'
likelihood <- function(tree,seqs,method='pruning',workenv=NULL) {
  tree=unroot(tree)
  if (is.null(workenv)) workenv=list()#use a list, not an environment

    if (method=='pruning') {
      if (is.null(workenv$pml)) {
        workenv$pml=phangorn::pml(tree,seqs,model='JC')
      } else {
        workenv$pml=stats::update(workenv$pml,tree=tree)
      }
      return(workenv$pml$logLik)
    }

  if (method=='joint' || method=='marginal') {
    #Using joint (as in eg Pupko et al. 2000) or marginal ancestral state reconstruction
    tree$node.label=paste0('N',1:tree$Nnode)
    if (is.null(workenv$pml)) {
      workenv$pml=phangorn::pml(tree,seqs,model='JC')
    } else {
      workenv$pml$tree=tree
    }
    if (method=='joint') {
      anc=phangorn::anc_pml(workenv$pml)$state
    } else {
      anc=phangorn::ancestral.pml(workenv$pml,type = "ml",return='ancestral')$state
    }
    r=likelihood_anc(tree,seqs,anc,(if (is.list(workenv)) NULL else workenv))
    return(r)
  }

  if (method=='parsimony') {
    #testing with parsimony ancestral reconstruction
    tree$node.label=paste0('N',1:tree$Nnode)
    anc=phangorn::ancestral.pars(tree,seqs,return='ancestral')
    anc=anc[tree$node.label]
    r=likelihood_anc(tree,seqs,anc)
    return(r)
  }
}

#' EXPERIMENTAL MLE inference of the tree given the sequences
#' @param seqs Sequences in phangorn::phyDat format
#' @returns Phylogenetic tree estimated
#'
experimental_mle <- function(seqs) {
  random_start=rtree(length(seqs),tip.label=names(seqs))
  random_start$edge.length=random_start$edge.length*0.01
  o=capture.output(m<-phangorn::pml_bb(seqs,model='JC',start=unroot(random_start)))
  return(m$tree)
}

#' EXPERIMENTAL Bayesian inference of the tree given the sequences
#' @param seqs Sequences in phangorn::phyDat format
#' @param niter Number of iterations
#' @param thin Thinin value
#' @param method Which method to use: pruning, joint, marginal or parsimony
#' @param showProgress Whether of not to show progress bar
#' @returns Result of Bayesian inference
#'
experimental_mcmc <- function(seqs,niter=10000,thin=100,method='pruning',showProgress=T) {
  proposal_probs  = c(NNI = 0.4, SPR = 0.4, branch = 0.2)
  tree=rtree(length(seqs),tip.label=names(seqs))
  tree$edge.length=tree$edge.length*0.01
  workenv=new.env()
  current_loglik=likelihood(tree,seqs,method=method,workenv=workenv)

  samples=vector("list",niter/thin)
  lik_trace=numeric(niter)
  acceptance=0
  if (showProgress) pb <- utils::txtProgressBar(min=0,max=niter,style = 3)
  for (i in 1:niter) {
    move_type=sample(names(proposal_probs), 1, prob = proposal_probs)
    tree_prop=tree
    if (move_type == "NNI") tree_prop=phangorn::rNNI(tree)
    if (move_type == "SPR") tree_prop=phangorn::rSPR(tree)
    if (move_type == "branch") tree_prop$edge.length=tree$edge.length * exp(rnorm(length(tree$edge.length), 0, 0.1))
    prop_loglik=likelihood(tree_prop,seqs,method=method,workenv=workenv)

    if (is.binary(tree_prop) && log(runif(1)) < prop_loglik - current_loglik) {
      tree=tree_prop
      current_loglik=prop_loglik
      acceptance=acceptance+1
    }

    lik_trace[i]=current_loglik
    if (i %% thin == 0) {
      if (showProgress) utils::setTxtProgressBar(pb, i)
      samples[[i/thin]]=tree
    }
  }
  if (showProgress) close(pb)

  return(list(samples=samples,lik_trace=lik_trace,acceptance=acceptance/niter))
}

#' Likelihood function with recombination given ancestral sequences
#' @param tree Phylogeny in ape::phylo format
#' @param seqs Sequences in phangorn::phyDat format
#' @param anc Ancestral sequences in phangorn::phyDat format
#' @param rho Rho parameter
#' @param delta Delta parameter
#' @param nu Nu parameter
#' @returns Log-likelihood
#'
likelihood_anc_rec <- function(tree,seqs,anc,rho,delta,nu) {
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

resolve_tip_uncertainties<-function(tree,seqs,anc)
{
  n=Ntip(tree)
  for (i in 1:nrow(tree$edge)) {
    if (tree$edge[i,2]>n) next
    bottom=seqs[[tree$tip.label[tree$edge[i,2]]]]
    top=anc[[tree$node.label[tree$edge[i,1]-n]]]
    inds=which(bottom==16)#missing data
    bottom[inds]=top[inds]
    seqs[[tree$tip.label[tree$edge[i,2]]]]=bottom
  }
  return(seqs)
}
