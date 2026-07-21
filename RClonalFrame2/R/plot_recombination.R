#' Plotting tree and recombination events of a given iteration from cf2out objects
#' @param x Output from ClonalFrame2
#' @param iteration Index of the iteration to plot (default is the last iteration)
#' @param gap Gap between sites in output (default is 100)
#' @param ... Additional parameters are passed on
#' @return Plot of ClonalFrame2 result
#'
plot_recombination = function(x,iteration,gap,...) {
  stopifnot(inherits(x, "cf2out"))
  if (missing(iteration)) iteration=nrow(x$iterations)
  if (missing(gap)) gap=100
  r=decode_recombination(x,iteration,gap)

}

