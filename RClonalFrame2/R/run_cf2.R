#' Run ClonalFrame2
#' @param inputfile Name of input FASTA file
#' @param outputfile Name of output CF2 file
#' @param iter Total number of iterations
#' @param thin Thinning interval
#' @param seed Seed for random number generator
#' @param ncpu Number of cores to use
#' @param profile Whether to activate the profiler
#' @param verbose Whether to activate verbose mode
#' @param simu_pars Vector containing the seven simulation parameters
#' @param starting_pars Vector containing the three starting parameters
#' @param update_pars Vector containing the five per-iteration update numbers
#' @param prior_pars Vector containing the eight prior parameters
#' @param path Path to ClonalFrame2 executable
#' @returns Object of class cf2out
#' @export
#'
run_cf2 <- function(inputfile,outputfile,iter,thin,seed,ncpu,profile,verbose,simu_pars,starting_pars,update_pars,prior_pars,path) {
  if (missing(inputfile) && missing(simu_pars)) stop('Need to provide either an input file or simulation parameters')
  if (!missing(inputfile) && !missing(simu_pars)) stop('Cannot use both an input file and simulation parameters')
  if (!missing(inputfile) && !file.exists(inputfile)) stop('Input file is not found')
  if (missing(path)) path='ClonalFrame2'
  path=Sys.which(path)
  if (nchar(path)==0) stop('Could not find ClonalFrame2 executable')
  args=c()
  if (!missing(inputfile)) args=c(args,'-i',shQuote(path.expand(inputfile)))
  if (missing(outputfile)) outputfile=tempfile()
  args=c(args,'-o',shQuote(outputfile))
  if (!missing(iter)) args=c(args,'-x',iter)
  if (!missing(thin)) args=c(args,'-y',thin)
  if (!missing(seed)) args=c(args,'-s',seed) else args=c(args,'-s',round(runif(1,min=0,max=1e9)))
  if (!missing(ncpu)) args=c(args,'-c',ncpu)
  if (!missing(profile) && profile) args=c(args,'-p')
  if (!missing(verbose) && verbose) args=c(args,'-v')
  if (!missing(simu_pars)) args=c(args,'-S',sprintf('%d,%d,%f,%f,%f,%f,%f',simu_pars[1],simu_pars[2],simu_pars[3],simu_pars[4],simu_pars[5],simu_pars[6],simu_pars[7]))
  if (!missing(starting_pars)) args=c(args,'-I',sprintf('%f,%f,%f',starting_pars[1],starting_pars[2],starting_pars[3]))
  if (!missing(update_pars)) args=c(args,'-U',sprintf('%d,%d,%d,%d,%d',update_pars[1],update_pars[2],update_pars[3],update_pars[4],update_pars[5]))
  if (!missing(prior_pars)) args=c(args,'-P',sprintf('%f,%f,%f,%f,%f,%f,%f,%f',prior_pars[1],prior_pars[2],prior_pars[3],prior_pars[4],prior_pars[5],prior_pars[6],prior_pars[7],prior_pars[8]))
  message('Running ClonalFrame2...')
  o=suppressWarnings(file.remove(outputfile))
  status=system2(path,args)
  if (status!=0) stop(sprintf('ClonalFrame2 exited with failed status %d',status))
  if (!file.exists(outputfile)) stop('Output file was not generated')
  message('Reading ClonalFrame2 output...')
  load_cf2(outputfile)
}
