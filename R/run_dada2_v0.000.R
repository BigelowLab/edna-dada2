#' Process using dada2 with automated <something>
#' 
#' 
#' 
library(dada2)
library(yaml)
library(dplyr)
library(readr)
library(ShortRead)  
library(Biostrings)
library(configr)
library(futile.logger)


#' Given a path - make it if it doesn't exist
#' 
#' We don't use R's \code{dir.create} in case we a reaching across volumes which 
#' has caused headaches at time.
#' 
#' @param path character, the path to check and/or create
#' @param logical, TRUE if the path exists
make_path <- function(path){
  ok <- dir.exists(path[1])
  if (!ok){
    ok <- system2("mkdir", args = path[1]) == 0
  }
  ok
}

#' Count the number of matches of the primer in the subject
#' 
#' @param primer character, pattern string
#' @param fn XString subject to match
#' @return numeric, count of hits
count_primer_hits <- function(primer, fn) {
  nhits <- Biostrings::vcountPattern(primer, 
                                     ShortRead::sread(ShortRead::readFastq(fn)), 
                                     fixed = FALSE)
  return(sum(nhits > 0))
}


#' Given a primer, compute all of the possible orientations
#' 
#' @param primer input sequence, see \code{\link[Biostrings]{DNAString}}
#' @return name character vector
all_orients <- function(primer = "TTGAAAA-CTC-N")  {                                                                                                                                 
  dna <- Biostrings::DNAString(primer)                                                                                                           
  orients <- c(Forward = dna, 
             Complement =  complement(dna), 
             Reverse = reverse(dna), 
             RevComp = reverseComplement(dna)) 
  return(sapply(orients, toString))                                                                                                 
} 

#' Check the configuration to make sure it is complete, borrowing from defaults as needed.
#'
#'
#' @param x list, configuration list
#' @param default, list default configuration
check_configuration <- function(
	x = list(foo = list(bar = 7, biz = "h"), 
	         dog = list(drool = TRUE)), 
	default = list(foo = list(bar = 9, biz = "p"), 
								 fish = list(pickled = "yum", breakfast = TRUE),
	               dog = list(drool = TRUE, quantity = "a lot"),
	               cat = list(hairball = TRUE))){

	nm_d <- names(default)
	ix <- nm_d %in% names(x)
	if (any(!ix))	x[nm_d[!ix]] <- default[nm_d[!ix]]
	
	for (nm in nm_d){
		nms <- names(default[[nm]])
		ix <- nms %in% names(x[[nm]])
		if (any(!ix)){
			wix <- which(ix)
			for (i in wix) x[[nm]][[nms[i]]] <- default[[nm]][[nms[i]]]
		}
	}
	x
}


#' Default configuration - softwired for now
#'
#' @return named and nested configuration list
default_configuration <- function(){

	list(
	  email = "btupper@bigelow.org", 
    input_path = ".", 
    output_path = ".", 
    dada2 = list(
      maxN = 0, 
      multithread = 32, 
      compress = FALSE), 
    cutadapt = list(
    	app = "/mnt/modules/bin/dada2/1.12/bin/cutadapt", 
      `--minimum-length` = 25, 
      `-n` = 2), 
    primer = list(
        FWD = "CYGCGGTAATTCCAGCTC", 
        REV = "AYGGTATCTRATCRTCTTYG"))
}

#' Get the configuration use default if needed
#' 
#' @param x configuration filename
#' @param default list describing the default
#' @return list config values
get_configuration <- function( x =  commandArgs(trailingOnly = TRUE),
  default = default_configuration()){
  
  if (length(x) == 0) {
    cfg <- default_configuration()
  } else {
    cfg <- try(configr::read.config(x[1]))
    if (inherits(cfg, 'try-error')){
      print(cfg)
      stop("failed to read config file:", x[1])
    }
    # TODO
    cfg <- check_configuration(cfg, default)
  }
  cfg
}


#' List fastq files and separate into forward and reverse reads
#' 
#' @param path character, the input path
#' @param pattern_forward file pattern 
#' @param pattern_reverse file pattern
#' @return named list of sorted foreward and reverse fastq filenames
list_fastq <- function(path,
                       pattern_forward = "^.*_R1_001\\.fastq",
                       pattern_reverse = "^.*_R2_001\\.fastq"){
  
  list(
    forward = sort(list.files(path, pattern = pattern_forward, full.names = TRUE)),
    reverse = sort(list.files(path, pattern = pattern_reverse, full.names = TRUE)) )
}

#' Filter and trim
#' 
#' @param fq list of forward and reverse fastq files
#' @param subdirectory the directory where the new fastqs are stored
#' @param maxN numeric see \code{\link[dada2]{filterAndTrim}}
#' @param multithread numeric see \code{\link[dada2]{filterAndTrim}}
#' @param compress logical see \code{\link[dada2]{filterAndTrim}}
#' @param ... other arguments for \code{\link[dada2]{filterAndTrim}}
#' @return integer matrix see \code{\link[dada2]{filterAndTrim}}
filter_and_trim <- function(fq, 
                            subdirectory = 'filtN',
                            maxN = 0, 
                            multithread = 32, 
                            compress = FALSE,
                            ...){
  
  ffilt <- file.path(dirname(fq$forward), subdirectory, basename(fq$forward))
  rfilt <- file.path(dirname(fq$reverse), subdirectory, basename(fq$reverse))
  dada2::filterAndTrim(fq$forward, ffilt, 
                               rev = fq$reverse, filt.rev = rfilt,
                               maxN = maxN, 
                               multithread = multithread, 
                               compress = compress,
                               ...)
}



# main processing step
main <- function(){

CFG <- get_configuration(x = "/home/btupper/edna/edna-dada2/config/run_dada2_v0.000.yml")

if (nchar(Sys.which(CFG$cutadapt$app)) == 0){
  stop("cutadapt application not found:", CFG$cutadapt$app)
}

if (!dir.exists(CFG$input_path)) stop("input path not found:", CFG$input_path)
if (!make_path(CFG$output_path)) stop("output path not created:", CFG$output_path)

fq_files <- list_fastq(CFG$input_path)

if (length(fq_files[[1]]) <= 0)
    stop("fastq files not found:", CFG$input_path )
if (length(fq_files[[1]]) != length(fq_files[[2]]))
    stop(sprintf("unequal number of fastq files: %i forward and %i reverse", 
                 length(fq_files$forward), length(fq_files$reverse)))


FWD.orients <- all_orients(CFG$primer$FWD) 
REV.orients <- all_orients(CFG$primer$REV) 

mat <- filter_and_trim(fq_files, 
                       subdirectory = "filtN",
                       maxN = CFG$dada2$maxN, 
                       multithread = CFG$dada2$multithread, 
                       compress = CFG$dada2$compress)

}




