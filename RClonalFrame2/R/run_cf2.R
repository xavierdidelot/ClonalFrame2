#' Run ClonalFrame2
#' @param inputfile Name of input FASTA file
#' @param outputfile Name of output CF2 file
#' @param args Optional arguments to use in CF2 run
#' @param path Path to ClonalFrame2 executable
#' @returns Object of class cf2out
#' @export
#'
run_cf2 <- function(inputfile,outputfile,args,path) {
  if (missing(inputfile)) stop('Need to provide inputfile')
  if (!file.exists(inputfile)) stop('Input file is not found')
  if (missing(outputfile)) outputfile='/tmp/out.json'
  if (missing(args)) args=''
  if (missing(path)) path='/Users/u1775021/github/ClonalFrame2/src'
  message('Compiling if needed...')
  system(sprintf('make -C %s',path))
  message('Running ClonalFrame2...')
  system(sprintf('%s/ClonalFrame2 -i %s -o %s %s',path,inputfile,outputfile,args))
  message('Reading ClonalFrame2 output...')
  load_cf2(outputfile)
}
