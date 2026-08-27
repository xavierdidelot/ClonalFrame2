#' Run ClonalFrame2
#' @param inputfile Name of input FASTA file
#' @param outputfile Name of output CF2 file
#' @param iter Total number of iterations
#' @param thin Thining interval
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
  args=''
  if (!missing(inputfile)) args=sprintf('%s -i %s ',args,inputfile)
  if (missing(outputfile)) outputfile=tempfile()
  args=sprintf('%s -o %s ',args,outputfile)
  if (!missing(iter)) args=sprintf('%s -x %d ',args,iter)
  if (!missing(thin)) args=sprintf('%s -y %d ',args,thin)
  if (!missing(seed)) args=sprintf('%s -s %d ',args,seed) else args=sprintf('%s -s %d ',args,round(runif(1,min=0,max=1e9)))
  if (!missing(ncpu)) args=sprintf('%s -c %d ',args,ncpu)
  if (!missing(profile) && profile) args=sprintf('%s -p ',args)
  if (!missing(verbose) && verbose) args=sprintf('%s -v ',args)
  if (!missing(simu_pars)) args=sprintf('%s -S %d,%d,%f,%f,%f,%f,%f ',args,simu_pars[1],simu_pars[2],simu_pars[3],simu_pars[4],simu_pars[5],simu_pars[6],simu_pars[7])
  if (!missing(starting_pars)) args=sprintf('%s -I %f,%f,%f ',args,starting_pars[1],starting_pars[2],starting_pars[3])
  if (!missing(update_pars)) args=sprintf('%s -U %d,%d,%d,%d,%d ',args,update_pars[1],update_pars[2],update_pars[3],update_pars[4],update_pars[5])
  if (!missing(prior_pars)) args=sprintf('%s -P %f,%f,%f,%f,%f,%f,%f,%f ',args,prior_pars[1],prior_pars[2],prior_pars[3],prior_pars[4],prior_pars[5],prior_pars[6],prior_pars[7],prior_pars[8])
  message('Running ClonalFrame2...')
  o=suppressWarnings(file.remove(outputfile))
  system(sprintf('MallocNanoZone=;%s %s',path,args))
  if (!file.exists(outputfile)) stop('Output file was not generated')
  message('Reading ClonalFrame2 output...')
  load_cf2(outputfile)
}
