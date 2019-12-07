#' Process using dada2 with automated <something>
#' 
#' 
#' 
library(dada2)
library(yaml)
library(readr)
library(ShortRead)  
library(Biostrings)
library(config)
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
count_primer_its <- function(primer, fn) {
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

#' Get the configuration use default if needed
#' 
#' @param default list describing the default
#' @return list config values
get_configuration <- function(
  default = list(email = "btupper@bigelow.org", input+path = ".", output_path = "./output", 
                 dada2 = list(maxN = 0, multithread = 32, compress = FALSE), 
                 cutadapt = list(app = "/mnt/modules/bin/dada2/1.12/bin/cutadapt", 
                                 `--minimum-length` = 25, `-n` = 2))){
  
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0) {
    cfg <- default
  } else {
    cfg <- try(configr::read.config(args[1]))
    if (inherits(cfg, 'try-error'))
      print(cfg)
      stop("failed to read config file:", args[1])
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
                            subdiretcory = 'filtN',
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

CFG <- get_configuration()

if (nchar(Sys.which(CFG$cutadapt$app)) == 0){
  stop("cutadapt application not found:", CFG$cutadapt$app)
}

if (!file.exists(CFG$input_path)) stop("input path not found:", CFG$input_path)
if (!make_path(CFG$output_path)) stop("output path not created:", CFG$output_path)

fq <- list_fastq(CFG$input_path)
if (!all.equal(lengths(fq)))
    stop(sprintf("unequal number of fastq files: %i forward and %i reverse", 
                 length(fq$forward), length(fq$reverse)))

mat <- filter_and_trim(fq, 
                       subdirectory = "filtN",
                       maxN = CFG$dada2$maxN, 
                       multithread = CFG$dada2$multithread, 
                       compress = CFG$dada2$compress)





