#' Process using dada2 with automated <something>
#' 
#' 
#' 

library(dplyr)
library(readr)
library(configr)
library(futile.logger)
library(patchwork)

library(ShortRead)  
library(Biostrings)
library(dada2)



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
primer_hits <- function(primer, fn) {
  nhits <- Biostrings::vcountPattern(primer, 
                                     ShortRead::sread(ShortRead::readFastq(fn)), 
                                     fixed = FALSE)
  return(sum(nhits > 0))
}



#' Compute primer counts
#'
#' @param FWD character, forward primer
#' @param REV character, reverse primer
#' @param fn_FWD XString forward subject to match
#' @param fn_REV XString, reverse subject to match
#' @param form character, desired output format - either 'matrix' with row names or "table" (tibble)
#' @return table or matrix of counts
primer_counts <- function(FWD, REV, fn_FWD, fn_REV,
  form = c("matrix", "table")[2]){
  
  x <- rbind(
    FWD.ForwardReads = sapply(FWD, primer_hits, fn = fn_FWD[[1]]),
    FWD.ReverseReads = sapply(FWD, primer_hits, fn = fn_REV[[1]]),
    REV.ForwardReads = sapply(REV, primer_hits, fn = fn_FWD[[1]]),
    REV.ReverseReads = sapply(REV, primer_hits, fn = fn_REV[[1]]))
  
  if (tolower(form[1]) == 'table'){
    x <- dplyr::as_tibble(x, rownames = "name")
  }
  
  x
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
  if (any(!ix))  x[nm_d[!ix]] <- default[nm_d[!ix]]
  
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
      more_args = "--minimum-length = 25 -n = 2"), 
    primer = list(
        FWD = "CYGCGGTAATTCCAGCTC", 
        REV = "AYGGTATCTRATCRTCTTYG"),
    taxa_level = c("Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species"))
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


#' Strip the extension(s) off of a filename
#'
#' Note if ext is ".fastq" then ".fastq.gz" and ".fastq.tar.gz" will also be stripped
#''
#' @param filename character one or more filenames
#' @param ext character, one or more extension patterns
#' @return filename with extension stripped
strip_extension <- function(filename = c("BR2_2016_S216_L001_R2_001.fastq", "foobar.fastq.gz", "fuzzbaz.txt"),
  ext = ".fastq"){
  
  ix <- gregexpr(ext, filename, fixed = TRUE)
  sapply(seq_along(ix), 
    function(i){
      if (ix[[i]] != -1) {
        s <- substring(filename[i], 1, ix[[i]]-1)
      } else {
        s <- filename[i]
      }
      s
    })
  
  
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
#' @param form desired output format - either matrix with rown names of table (tibble)
#' @param ... other arguments for \code{\link[dada2]{filterAndTrim}}
#' @return integer matrix see \code{\link[dada2]{filterAndTrim}}
filter_and_trim <- function(fq, 
                            subdirectory = 'filtN',
                            form = c("matrix", "table")[2],
                            ...){
  
  ffilt <- file.path(dirname(fq$forward), subdirectory, basename(fq$forward))
  rfilt <- file.path(dirname(fq$reverse), subdirectory, basename(fq$reverse))
  x <- dada2::filterAndTrim(fq$forward, 
  													ffilt, 
                            rev = fq$reverse, 
                            filt.rev = rfilt,
                            ...)
  if (tolower(form[1]) == 'table') {
    n <- rownames(x)
    x <- dplyr::as_tibble(x, rownames = "name")
  }
  x
}


#' Run cutadapt
#'
#' @param cut_files two element list for forward and reverse cutadapt results to be created
#' @param filt_files two element list of forward and reverse filtered files
#' @param CFG list of configuration
#' @param save_output logical, if TRUE try to capture the output of cutadapt to text files,
#'   one per cutadapt file generated, in the cutadapt outpuyt directory
#' @param save_graphics logical, if TRUE try to capture quality plots form the resulting cut_files
#' @return numeric codes, one per cut_file pairing as per \code{system2) where 0 means success'
run_cutadapt <- function(
  cut_files,
  filt_files,
  CFG,
  save_output = TRUE,
  save_graphics = FALSE){
  
  FWD.RC <- dada2:::rc(CFG$primer$FWD)
  REV.RC <- dada2:::rc(CFG$primer$REV)
  R1.flags <- paste("-g", CFG$primer$FWD, "-a", REV.RC)
  R2.flags <- paste("-G", CFG$primer$REV, "-A", FWD.RC)

  OK <- sapply(seq_along(cut_files$forward),
    function(i){
      if (save_output){
        ofile <- paste0(strip_extension(cut_files$forward[1]),".cutadapt_output.txt")
      } else {s
        ofile <- ""
      }
      system2(CFG$cutadapt$app, 
        args = c(
          R1.flags, 
          R2.flags, 
          CFG$more_args, 
           "-o", cut_files$forward[i], 
           "-p", cut_files$reverse[i],
           filt_files$forward[i], 
           filt_files$reverse[i]),
         stdout = ofile)
    })
  if (all(OK == 0) & save_graphics){
  	ix <- seq_len(max(length(cut_files$forward), 2))
  	ofile <- paste0(strip_extension(cut_files$forward[1]),".cutadapt_quality.pdf")
    pdf(ofile)
  	try(dada2::plotQualityProfile(cut_files$forward[ix]) +  dada2::plotQualityProfile(cut_files$reverse[ix]))
  	dev.off()
  }
  OK
}




# main processing step
main <- function(x = "/home/btupper/edna/edna-dada2/config/run_dada2_v0.000.yml"){

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

filtN_r <- filter_and_trim(fq_files, 
                       subdirectory = CFG$dada2_trimAndFilter_filtN$subdirectory,
                       maxN = CFG$dada2_trimAndFilter_filtN$maxN, 
                       multithread = CFG$dada2_trimAndFilter_filtN$multithread, 
                       compress = CFG$dada2_trimAndFilter_filtN$compress)

filtN_files <- list_fastq(file.path(CFG$input_path, "filtN"))

pcounts <- primer_counts(FWD.orients, REV.orients, filtN_files$forward, filtN_files$reverse)

cut_path <- file.path(CFG$input_path, "cutadapt")
if (!make_path(cut_path)) stop("cut_path not created:", cut_path)

cut_files <- lapply(filt_files,
  function(f){
    file.path(cut_path, basename(f))
  })

# turn of graphics until we find issue with printing reverse
cut_ok <- run_cutadapt(cut_files, filtN_files, CFG, save_output = TRUE,, save_graphcs = FALSE)
cut_pcounts <- primer_counts(FWD.orients, REV.orients, cut_files$forward, cut_files$reverse)

filtered_r <- filter_and_trim(cut_files, 
                       subdirectory = CFG$dada2_trimAndFilter_filtered$subdirectory,
                       trunQ = CFG$dada2_trimAndFilter_filtered$trunQ,
                       trunLen = CFG$dada2_trimAndFilter_filtered$trunLen,
                       maxEE = CFG$dada2_trimAndFilter_filtered$maxEE,
                       rm.phix= CFG$dada2_trimAndFilter_filtered$rm.phix,
                       maxN = CFG$dada2_trimAndFilter_filtered$maxN, 
                       multithread = CFG$dada2_trimAndFilter_filtered$multithread, 
                       compress = CFG$dada2_trimAndFilter_filtered$compress)





} #main


