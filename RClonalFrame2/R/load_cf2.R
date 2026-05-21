#' Load a CF2 output file
#' @param filename Name of output CF2 file
#' @returns Object of class cf2out
#' @export
#'
load_cf2 <- function(filename) {
  x=fromJSON(filename)
  class(x)<-"cf2out"
  return(x)
}

#' Print function for cf2out objects
#' @param x Output from ClonalFrame2
#' @param ... Passed on
#' @return Print out details of ClonalFrame2 results
#' @export
print.cf2out <- function(x, ...) {
  stopifnot(inherits(x, "cf2out"))
  cat(sprintf('ClonalFrame2 output with %d iterations\n',nrow(x$iterations)))
  invisible(x)
}

#' Summary function for cf2out objects
#' @param x Output from ClonalFrame2
#' @param ... Passed on
#' @return Print out details of ClonalFrame2 results
#' @export
summary.cf2out <- function(object, ...) {
  stopifnot(inherits(object, "cf2out"))
  cat(sprintf('ClonalFrame2 output with %d iterations\n',nrow(object$iterations)))
  invisible(object)
}

#' Plotting method for cf2out objects
#' @param x Output from ClonalFrame2
#' @param type Type of plot to do
#' @param ... Additional parameters are passed on
#' @return Plot of ClonalFrame2 results
#' @export
plot.cf2out = function(x, type='trace', ...) {
  old.par=par(no.readonly = T)

  if (type=='trace') {
    par(mfrow=c(2,2))
    plot(x$iterations[,'likelihood'],main='Likelihood',type='l',xlab='Sampled iterations',ylab='')
    plot(x$iterations[,'rho'],main='Rho',type='l',xlab='Sampled iterations',ylab='')
    plot(x$iterations[,'delta'],main='Delta',type='l',xlab='Sampled iterations',ylab='')
    plot(x$iterations[,'nu'],main='Nu',type='l',xlab='Sampled iterations',ylab='')
  }

  if (type=='trace_acceptance') {
    par(mfrow=c(2,2))
    plot(x$iterations[,'paccepted_tree'],main='Tree',type='l',xlab='Sampled iterations',ylab='')
    plot(x$iterations[,'paccepted_rho'],main='Rho',type='l',xlab='Sampled iterations',ylab='')
    plot(x$iterations[,'paccepted_delta'],main='Delta',type='l',xlab='Sampled iterations',ylab='')
    plot(x$iterations[,'paccepted_nu'],main='Nu',type='l',xlab='Sampled iterations',ylab='')
  }

  if (type=='trace_trees') {
    showing=4
    par(mfrow=c(ceiling(sqrt(showing)),ceiling(sqrt(showing))))
    inds=round(seq(1,nrow(x$iterations),length.out=showing))
    for (i in 1:showing) {
      t=read.tree(text=x$iterations[inds[i],'tree'])
      plot(t)
      axisPhylo(1)
    }
  }

  par(old.par)
}
