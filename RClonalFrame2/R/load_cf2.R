#' Load a CF2 output file
#' @param filename Name of output CF2 file
#' @returns Object of class cf2out
#' @export
#'
load_cf2 <- function(filename) {
  x=fromJSON(filename)
  x$pcburn=x$mcmc_parameters$burn/x$mcmc_parameters$nbit
  class(x)<-"cf2out"
  return(x)
}

#' Set the percentage of burnin
#' @param x Output from ClonalFrame2
#' @param val New value for percentage of burnin
#' @export
pcburn <- function(x, val) {
  stopifnot(inherits(x, "cf2out"))
  x$pcburn=val
  return(x)
}

#' Print function for cf2out objects
#' @param x Output from ClonalFrame2
#' @param ... Passed on
#' @return Print out details of ClonalFrame2 results
#' @export
print.cf2out <- function(x, ...) {
  stopifnot(inherits(x, "cf2out"))
  cat(sprintf('ClonalFrame2 output with %d sampled iterations and %d%% burnin\n',nrow(x$iterations),100*x$pcburn))
  iters=max(1,floor(nrow(x$iterations)*x$pcburn)):nrow(x$iterations)
  for (nam in c('rho','delta','nu')) {
    v=x$iterations[iters,nam]
    v=sort(v)
    vals=c(mean(v),v[pmax(1,floor(length(v)*c(0.025,0.975)))])
    cat(sprintf('%s=%.2e [%.2e;%.2e]\n',nam,vals[1],vals[2],vals[3]))
  }
  invisible(x)
}

#' Summary function for cf2out objects
#' @param object Output from ClonalFrame2
#' @param ... Passed on
#' @return Print out details of ClonalFrame2 results
#' @export
summary.cf2out <- function(object, ...) {
  stopifnot(inherits(object, "cf2out"))
  print(object,...)
  invisible(object)
}

#' Plotting method for cf2out objects
#' @param x Output from ClonalFrame2
#' @param type Type of plot to do
#' @param ... Additional parameters are passed on
#' @return Plot of ClonalFrame2 results
#' @export
plot.cf2out = function(x, type='pars', ...) {
  old.par=par(no.readonly = T)
  iters=max(1,floor(nrow(x$iterations)*x$pcburn)):nrow(x$iterations)

  if (type=='pars') {
    panel.hist <- function(x,truths,...)
    {
      usr <- par("usr")
      par(usr = c(usr[1:2], 0, 1.5) )
      h <- hist(x, breaks = seq(min(x),max(x),length.out=21),plot=F)
      breaks <- h$breaks; nB <- length(breaks)
      y <- h$counts; y <- y/max(y)
      rect(breaks[-nB], 0, breaks[-1], y, ...)
      pairs_ind<<-pairs_ind+1
      if (!is.null(truths)) lines(c(truths[pairs_ind],truths[pairs_ind]),c(0,1.5),col='red')
    }
    varnames=c("rho", "delta", "nu")
    truth_values=NULL
    if (!is.null(x$simulation_parameters)) truth_values=x$simulation_parameters[varnames]
    pairs_ind<<-0
    pairs(x$iterations[iters,varnames],diag.panel=function (x,...) panel.hist(x,truths=truth_values,...))
  }

  if (type=='trace') {
    par(mfrow=c(3,3))
    plot(iters,x$iterations[iters,'likelihood'],main='Likelihood',type='l',xlab='Sampled iterations',ylab='')
    plot(iters,x$iterations[iters,'prior_tree'],main='Tree prior',type='l',xlab='Sampled iterations',ylab='')
    plot(iters,x$iterations[iters,'prior_rho']+x$iterations[iters,'prior_delta']+x$iterations[iters,'prior_nu'],main='Params prior',type='l',xlab='Sampled iterations',ylab='')
    plot(iters,x$iterations[iters,'likelihood']+x$iterations[iters,'prior_tree']+x$iterations[iters,'prior_rho']+x$iterations[iters,'prior_delta']+x$iterations[iters,'prior_nu'],main='Posterior',type='l',xlab='Sampled iterations',ylab='')
    plot(iters,x$iterations[iters,'sum_bra_len'],main='Sum branch lengths',type='l',xlab='Sampled iterations',ylab='')
    plot(iters,x$iterations[iters,'rho'],main='Rho',type='l',xlab='Sampled iterations',ylab='')
    plot(iters,x$iterations[iters,'delta'],main='Delta',type='l',xlab='Sampled iterations',ylab='')
    plot(iters,x$iterations[iters,'nu'],main='Nu',type='l',xlab='Sampled iterations',ylab='')
  }

  if (type=='trace_acceptance') {
    iters=setdiff(iters,1)#first iteration always has zero paccept since nothing proposed yet
    par(mfrow=c(2,3))
    plot(x$iterations[iters,'paccepted_rho'],main='Rho',type='l',xlab='Sampled iterations',ylab='')
    plot(x$iterations[iters,'paccepted_delta'],main='Delta',type='l',xlab='Sampled iterations',ylab='')
    plot(x$iterations[iters,'paccepted_nu'],main='Nu',type='l',xlab='Sampled iterations',ylab='')
    plot(x$iterations[iters,'paccepted_lengths'],main='Lengths',type='l',xlab='Sampled iterations',ylab='')
    plot(x$iterations[iters,'paccepted_spr'],main='SPR',type='l',xlab='Sampled iterations',ylab='')
  }

  if (type=='trace_trees') {
    showing=4
    par(mfrow=c(ceiling(sqrt(showing)),ceiling(sqrt(showing))))
    inds=iters[round(seq(1,length(iters),length.out=showing))]
    for (i in 1:showing) {
      t=read.tree(text=x$iterations[inds[i],'tree'])
      plot(t)
      axisPhylo(1)
    }
  }

  if (type=='simulation') {
    plot_recombination(x,...)
  }

  if (type=='simulation_with_subs') {
    plot_recombination(x,showSubs=T,...)
  }

  par(old.par)
}
