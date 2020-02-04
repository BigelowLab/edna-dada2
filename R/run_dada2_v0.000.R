#' Process using dada2 with automated <something, threshold selection?>
#' 
#' Usage: from shell...
#' $ Rscript [opts] script_file config_file
#' 
#' @param [opts] any of the R invocation options described here
#'         https://cran.r-project.org/doc/manuals/r-release/R-intro.html#Invoking-R
#'         Suggested options are --no-save --no-restore --no-site-file --no-environ
#'         which will prevent headaches when R is launched and some damn fool saved
#'         previous R session info (Ben!)
#' @param script_file the full path and filename to the script.  For example
#'        /mnt/edna/code/run_dada2_v0.000.Rscript (please avoid spaces in paths!)  See an
#'        example here ... https://github.com/BigelowLab/edna-dada2/blob/master/R/run_dada2_v0.000.R
#' @param config_file the full path and filename for the simple-text YAML file. 
#'        The config file must have as a minimum entries for email, input_path and
#'        output_path. Any other entries should have default values embedded in the 
#'        script.  See and example here... 
#'        https://github.com/BigelowLab/edna-dada2/blob/master/config/run_dada2_v0.000.yml
#'
#' [NOTE]
#' Input requirements are that the input_path defined in the config file MUST have
#' complementary forward and reverse FASTQ files with the following naming patterns.
#'  	pattern_forward = "^.*_R1_001\\.fastq"
#'    pattern_reverse = "^.*_R2_001\\.fastq"
#' This seems quite strict to me, but I don't know much about fastq - is it compressed?  If not
#' should we allow for patterns like "^.*_R1_001\\.fastq.gz"?  Whatever the answer, we could
#' place the allowable patterns in the config file, but I think it would be better to simply
#' make the list_fastq() function more robust.
 
library(dplyr)
library(readr)
library(configr)
library(futile.logger)
library(patchwork)

library(ShortRead)  
library(Biostrings)
library(dada2)

#' Count the number of CPUs
#'
#' THis is a wrapper around \code{\link[parallel]{detectCores}}
#'
#' @return integer
count_cores <- function(){
  parallel::detectCores()
}

#' Given a path - make it if it doesn't exist
#' 
#' We don't use R's \code{dir.create} in case we a reaching across volumes which 
#' has caused headaches at time.
#' 
#' @param path character, the path to check and/or create
#' @return logical, TRUE if the path exists
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

#' Get the configuration use default if needed
#' 
#' @param x configuration filename
#' @param default list describing the default
#' @return list config values
get_configuration <- function( x =  commandArgs(trailingOnly = TRUE),
  default = configr::read.config("/home/btupper/edna/edna-dada2/config/run_dada2_v0.000.yml")){
  
  if (length(x) == 0) {
    cfg <- default
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
#' @param filelist list of forward and reverse fastq files
#' @param output_path character, the output path
#' @param maxN numeric see \code{\link[dada2]{filterAndTrim}}
#' @param multithread numeric see \code{\link[dada2]{filterAndTrim}}
#' @param compress logical see \code{\link[dada2]{filterAndTrim}}
#' @param ... other arguments for \code{\link[dada2]{filterAndTrim}}
#' @pram save_result logical, save CSV if TRUE
#' @return integer matrix as tibble see \code{\link[dada2]{filterAndTrim}}
filter_and_trim <- function(filelist, 
                            output_path = file.path(dirname(filelist$forward[1]),'filterAndTrim'),
                            save_results = FALSE,
                            ...){
  
  ffilt <- file.path(output_path, basename(filelist$forward))
  rfilt <- file.path(output_path, basename(filelist$reverse))
  x <- dada2::filterAndTrim(filelist$forward, 
                            ffilt, 
                            rev = filelist$reverse, 
                            filt.rev = rfilt,
                            ...)

  x <- dplyr::as_tibble(x, rownames = "name")
  if (save_results) {
    x <-  x %>%
    readr::write_csv(file.path(output_path, "filter_and_trim-results.csv"))
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


#' Run dada2::learnErrors on a set of fastq files
#'
#' @param filelist list of forward and reverse fastq files
#' @param ... arguments for \code{\link[dada2]{learnErrors}}
#' @param output_path character, the output path
#' @param save_output logical, if TRUE save the output to the specified output_path
#' @param save_graphics logical, if TRUE try to capture quality plots form the resulting cut_files
learn_errors <- function(filelist, ...,
  output_path = dirname(filelist$forward[1]),
  save_output = FALSE, 
  save_graphics = FALSE
  ){

  errs <- list(
      forward =  dada2::learnErrors(filelist$forward, ...),
      reverse =  dada2::learnErrors(filelist$forward, ...)
    )
    
  if (save_output){
    saveRDS(errs, file = file.path(output_path, "learn_errors.rds"))
  }
  
  if (save_graphics){
    pforward <- dada2::plotErrors(errs$forward, nominalQ=TRUE) + 
      ggplot2::ggtitle("Forward")
    preverse <- dada2::plotErrors(errs$reverse, nominalQ=TRUE) +
      ggplot2::ggtitle("Reverse")
    ofile <- file.path(output_path, "learn_errors.pdf")
    pdf(ofile, height = 7.5, width = 10.5)    
    try(
      print(pforward + preverse)
    )
    dev.off()
  }

  errs
}

#' Run dada2::dada
#'
#' @param filelist list of forward and reverse fastq files
#' @params errs list of forward and reverse outputs of learnErrors
#' @param ... arguments for \code{\link[dada2]{dada}}
run_dada <- function(filelist, errs, ...){
  

  filelist <- lapply(filelist, dada2::derepFastq)

  x <- list(
    forward = dada2::dada(filelist$forward, errs$forward, ...),
    reverse = dada2::dada(filelist$reverse, errs$reverse, ...)
  )
  x
}

#' Merge pairs ala dada2::mergePairs
#'
#' @param filelist list of forward and reverse fastq files
#' @param dada_r  list of dada2::dada results
merge_pairs <- function(filelist, dada_r, ...){
  filelist <- lapply(filelist, dada2::derepFastq)
  x <- dada2::mergePairs(
    dada_r$forward, 
    filelist$forward,
    dada_r$reverse, 
    filelist$reverse, 
    ...) 

  x
}

#' Count uniques ala \code{\link[dada2{getUniques}]}
#'
#' @param x object from which uniques-vector can be extracted
#' @param ... further arguments for \code{\link[dada2{getUniques}]}
#' @return 
count_uniques <- function(x, ...){
  sum(getUniques(x, ...)) 
}


#' main processing step - tries to operate as a pipeline returning 0 (success) or
#' failure ( > 0)
#' 
#' @param cfg character filename of the YAML configuration file
#' @return integer where 0 means success
main <- function(
  cfg = "/home/btupper/edna/edna-dada2/config/run_dada2_v0.000.yml"){


  RETURN <- 0
  CFG <- get_configuration(cfg)
  
  if (!dir.exists(CFG$input_path)){
    warning("input path not found:", CFG$input_path)
    return(RETURN + 1)
  }
  
  if (!make_path(CFG$output_path)){
    warning("output path not created:", CFG$output_path)
    return(RETURN + 1)
  }
  
  flog.threshold(toupper(CFG$verbose[1]))
  flog.appender(appender.tee(file.path(CFG$output_path, "log")) )
  flog.info("starting run: %s", cfg)
  

  if (is.numeric(CFG$dada2_dada_filtered$multithread)){
    MAX_CORES <- count_cores() - 1
	  CFG$dada2_dada_filtered$multithread <- pmax(CFG$dada2_dada_filtered$multithread, MAX_CORES)
  	flog.info("N cores: %i", MAX_CORES)
  }
  
  if (nchar(Sys.which(CFG$cutadapt$app)) == 0){
    flog.error("cutadapt application not found: %s", CFG$cutadapt$app)
    return(RETURN + 1)
  }
  
  flog.info("checking for input fastq files")
  fq_files <- list_fastq(CFG$input_path)
  sample.names <- sapply(strsplit(basename(fq_files$forward), "_"), `[`, 1)
  
  if (length(fq_files[[1]]) <= 0){
    flog.error("fastq files not found: %s", CFG$input_path)  
    return(RETURN + 1)
  }
  if (length(fq_files[[1]]) != length(fq_files[[2]])){
    flog.error("unequal number of fastq files: %i forward and %i reverse", 
               length(fq_files$forward), length(fq_files$reverse))
    return(RETURN + 1) 
  }
  
  flog.info("computing all orientations")
  FWD.orients <- all_orients(CFG$primer$FWD) 
  REV.orients <- all_orients(CFG$primer$REV) 
  
  ##### filterAndTrim filtN
  flog.info("filter and trim of input files")
  filtN_path <- file.path(CFG$output_path, CFG$dada2_filterAndTrim_filtN$name)
  if (!make_path(filtN_path)) {
    flog.error("filtN_path not created: %s", filtN_path)
    return(RETURN + 1)  
  }
  
  filtN_r <- filter_and_trim(fq_files, 
                         output_path = filtN_path,
                         maxN = CFG$dada2_filterAndTrim_filtN$maxN, 
                         multithread = CFG$dada2_filterAndTrim_filtN$multithread, 
                         compress = CFG$dada2_filterAndTrim_filtN$compress) %>%
    readr::write_csv(file.path(filtN_path, paste0(CFG$dada2_filterAndTrim_filtN$name, "-results.csv")))
  
  flog.info("primer counts")
  filtN_files <- list_fastq(filtN_path)
  pcounts <- primer_counts(FWD.orients, REV.orients, filtN_files$forward, filtN_files$reverse) %>%
      readr::write_csv(file.path(filtN_path, paste0(CFG$dada2_filterAndTrim_filtN$name, "-primer-counts.csv")))
  
  
  ###### cutadapt
  flog.info("cutadapt")
  cut_path <- file.path(CFG$output_path, "cutadapt")
  if (!make_path(cut_path)) {
    flog.error("cut_path not created:%", cut_path)
    return(RETURN + 1)
  }
  cut_files <- lapply(filtN_files,
    function(f){
      file.path(cut_path, basename(f))
    })
  
  # turn off graphics until we find issue with printing reverse
  cut_ok <- run_cutadapt(cut_files, filtN_files, CFG, save_output = TRUE, save_graphics = FALSE)
  if (all(cut_ok == 0)) {
    cut_pcounts <- primer_counts(FWD.orients, REV.orients, cut_files$forward, cut_files$reverse) %>%
      readr::write_csv(file.path(cut_path, paste0(basename(cut_path), "-primer-counts.csv")))
  } else {
    if (VERBOSE) flog.error("one or more cutadapt runs failed: %s", paste(cut_ok, collapse = ", "))
    return(RETURN + 1)
  }
  
  ##### filterAndTrim filtered
  flog.info("filter and trim of filtered")
  filtered_path <- file.path(CFG$output_path, CFG$dada2_filterAndTrim_filtered$name)
  if (!make_path(filtered_path)){ 
    flog.error("filtered_path not created: %s", filtered_path)
    return(RETURN + 1)
  }
  filtered_r <- filter_and_trim(cut_files, 
                         output_path = filtered_path,
                         truncQ = CFG$dada2_filterAndTrim_filtered$truncQ,
                         truncLen = CFG$dada2_filterAndTrim_filtered$truncLen,
                         maxEE = CFG$dada2_filterAndTrim_filtered$maxEE,
                         rm.phix= CFG$dada2_filterAndTrim_filtered$rm.phix,
                         maxN = CFG$dada2_filterAndTrim_filtered$maxN,  
                         compress = CFG$dada2_filterAndTrim_filtered$compress,
                         multithread = CFG$dada2_filterAndTrim_filtered$multithread) %>%
    readr::write_csv(file.path(filtered_path, paste0(CFG$dada2_filterAndTrim_filtN$name, "-results.csv")))
  
  filtered_files <- list_fastq(filtered_path)
  
  # learn errors
  flog.info("Learn errors")
  errs <- learn_errors(filtered_files,
    output_path = filtered_path,
    save_output = TRUE,
    save_graphics = TRUE)
  
  # run dada
  flog.info("run dada")
  dada_r <- run_dada(
    filtered_files, 
    errs, 
    multithread = CFG$dada2_dada_filtered$multithread)  # switch to CFG$dada2_dada_filtered
  
  # run merge pairs
  flog.info("merge pairs")
  mergers <- merge_pairs(filtered_files, dada_r, verbose = TRUE)
  saveRDS(mergers, file = file.path(filtered_path, "mergers.rds"))
  
  flog.info("make sequence table")
  seqtab <- dada2::makeSequenceTable(mergers) 
  tseqtab <- dplyr::as_tibble(t(seqtab)) %>%
    readr::write_csv(file.path(filtered_path, "seqtab.csv"))
  
  flog.info("remove Bimera Denovo")
  seqtab.nochim <- dada2::removeBimeraDenovo(seqtab, 
    method      = CFG$dada2_removeBimeraDenovo_seqtab$method, 
    multithread = CFG$dada2_removeBimeraDenovo_seqtab$multithread, 
    verbose     = CFG$dada2_removeBimeraDenovo_seqtab$verbose)
  tseqtab.nochim <- dplyr::as_tibble(t(seqtab.nochim)) %>%
    readr::write_csv(file.path(filtered_path, "seqtab-nochim.csv"))
  
  
  track <- dplyr::tibble(
                         name               = sample.names,
                         input              = filtered_r$reads.in, 
                         filtered           = filtered_r$reads.out,
                         denoised_forward   = sapply(dada_r$forward, count_uniques), 
                         denoised_reverse   = sapply(dada_r$reverse, count_uniques), 
                         merged             = sapply(mergers, count_uniques), 
                         nonchim            = rowSums(seqtab.nochim)) %>%
    readr::write_csv(file.path(filtered_path, "track.csv"))
  
  flog.info("assign taxonomy")
  taxa <- dada2::assignTaxonomy(seqtab.nochim, 
      refFasta          = CFG$dada2_assignTaxonomy_nochim$refFasta, 
      taxLevels         = CFG$dada2_assignTaxonomy_nochim$taxLevels, 
      minBoot           = CFG$dada2_assignTaxonomy_nochim$minBoot, 
      outputBootstraps  = CFG$dada2_assignTaxonomy_nochim$outputBootstraps, 
      verbose           = CFG$dada2_assignTaxonomy_nochim$verbose, 
      multithread       = CFG$dada2_assignTaxonomy_nochim$multithread) %>%
    dplyr::as_tibble() %>%
    read::write_csv(file.path(filtered_path, "taxa.csv"))

  return(RETURN)
} #main


# we only run is run as a script - not if interactive
if (!interactive()){
  cfgfile <- commandArgs(trailingOnly = TRUE)[1]
  ok <- main(cfgfile)
  quit(save = "no", status = ok)
}
