#' Run ClonalFrame2
#' @param inputfile Name of input FASTA file
#' @param outputfile Name of output CF2 file
#' @param burn Number of burn-in iterations
#' @param iter Total number of iterations
#' @param thin Thining interval
#' @param seed Seed for random number generator
#' @param ncpu Number of cores to use
#' @param profile Whether to activate the profiler
#' @param verbose Whether to activate verbose mode
#' @param simu_pars Vector containing the seven simulation parameters
#' @param update_pars Vector containing the four per-iteration update numbers
#' @param path Path to ClonalFrame2 executable
#' @returns Object of class cf2out
#' @export
#'
run_cf2 <- function(inputfile,outputfile,burn,iter,thin,seed,ncpu,profile,verbose,simu_pars,update_pars,path) {
  if (missing(inputfile) && missing(simu_pars)) stop('Need to provide either an input file or simulation parameters')
  if (!missing(inputfile) && !missing(simu_pars)) stop('Cannot use both an input file and simulation parameters')
  if (!missing(inputfile) && !file.exists(inputfile)) stop('Input file is not found')
  if (missing(path)) path='/Users/u1775021/github/ClonalFrame2/src'
  if (!file.exists(sprintf('%s/ClonalFrame2',path))) stop('Could not find ClonalFrame2 executable')
  args=''
  if (!missing(inputfile)) args=sprintf('%s -i %s ',args,inputfile)
  if (missing(outputfile)) outputfile='/tmp/out.json'
  args=sprintf('%s -o %s ',args,outputfile)
  if (!missing(burn)) args=sprintf('%s -x %d ',args,burn)
  if (!missing(iter)) args=sprintf('%s -y %d ',args,iter)
  if (!missing(thin)) args=sprintf('%s -z %d ',args,thin)
  if (!missing(seed)) args=sprintf('%s -s %d ',args,seed) else args=sprintf('%s -s %d ',args,round(runif(1,min=0,max=1e9)))
  if (!missing(ncpu)) args=sprintf('%s -c %d ',args,ncpu)
  if (!missing(profile) && profile) args=sprintf('%s -p ',args)
  if (!missing(verbose) && verbose) args=sprintf('%s -v ',args)
  if (!missing(simu_pars)) args=sprintf('%s -S %d,%d,%f,%f,%f,%f,%d ',args,simu_pars[1],simu_pars[2],simu_pars[3],simu_pars[4],simu_pars[5],simu_pars[6],simu_pars[7])
  if (!missing(update_pars)) args=sprintf('%s -U %d,%d,%d,%d ',args,update_pars[1],update_pars[2],update_pars[3],update_pars[4])
  message('Compiling if possible and needed...')
  if (file.exists(sprintf('%s/Makefile',path))) system(sprintf('make -C %s',path))
  message('Running ClonalFrame2...')
  system(sprintf('rm -rf %s',outputfile))
  system(sprintf('MallocNanoZone=;%s/ClonalFrame2 %s',path,args))
  if (!file.exists(outputfile)) stop('Output file was not generated')
  message('Reading ClonalFrame2 output...')
  load_cf2(outputfile)
}
