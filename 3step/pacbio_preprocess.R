library(dplyr)
library(readr)

library(ShortRead)  
library(Biostrings)
library(dada2)

library(charlier)
library(dadautils)


#' main processing step - tries to operate as a pipeline returning 0 (success) or
#' failure ( > 0)
#' 
#' @param CFG list, configuration from a config file
#' @return integer where 0 means success
main <- function(CFG){
  
  RETURN <- 0
  
  if (tolower(CFG$stage) != "preprocess"){
    warning("configuration stage should be 'preprocess'")
    return(RETURN + 1)
  }
  
  if (CFG$multithread[1] == "auto") CFG$multithread <- charlier::count_cores()
  
  if (!dir.exists(CFG$input_path)){
    warning("input path not found:", CFG$input_path)
    return(RETURN + 1)
  }
  
  if (!charlier::make_path(CFG$output_path)){
    warning("output path not created:", CFG$output_path)
    return(RETURN + 1)
  }
  
  # logging - see https://github.com/BigelowLab/charlier/wiki/Logging
  charlier::start_logger(filename = file.path(CFG$output_path, "log"))
  
  # audit - see https://github.com/BigelowLab/charlier/wiki/Auditing
  PBS_JOBID <- charlier::get_pbs_jobid(no_pbs_text = "not in PBS queue")
  charlier::audit(file.path(CFG$output_path, "audit.txt"), pbs_jobid = PBS_JOBID)
  
  # add baseline info into log, just because
  charlier::info("**** starting run ****")
  charlier::info("NCPUS: %s", as.character(CFG$multithread))
  charlier::info("System PID: %s", Sys.getpid())
  charlier::info("PBS_JOBID: %s", PBS_JOBID)
  charlier::info("VERSION: %s", CFG$version)
  charlier::info("INPUT PATH: %s", CFG$input_path)
  charlier::info("OUTPUT PATH: %s", CFG$output_path)
  
  charlier::info("checking for input fastq files")
  input_files <- auntie::list_filepairs(CFG$input_path, pattern_forward = "*.fastq.gz", verify=F) #%>%
  #dadautils::verify_filepairs()
  if (all( count_filepairs(input_files) == 0 )) {
    msg <- sprintf("check list_filepairs patterns, no files found in %s", CFG$input_path)
    charlier::error(msg)
    stop(msg)
  }
  
  # pacbio?
  norev <- (length(input_files) == 1) || (lengths(input_files)[[2]] == 0)
  
  sample.names <- dadautils::extract_sample_names(input_files, rule="basename")
  
  charlier::info("removing primers")
  remove_primers(input_files, primer.fwd=CFG$dada2_removePrimers$primer.fwd, 
                 primer.rev=dada2::rc(CFG$dada2_removePrimers$primer.rev), orient=CFG$dada2_removePrimers$orient)
  
  fq_files_no_primer <- auntie::list_filepairs(file.path(CFG$input_path, "no_primers"), 
                                                  verify=F,pattern_forward = "^.*.fastq.gz",pattern_reverse = "^.*.R2.fastq.gz")
  
  if("dada2_filterAndTrim" %in% names(CFG)){
    
    charlier::info("filter and trim of input files")
    
    filtN_path <- file.path(CFG$output_path, CFG$dada2_filterAndTrim$name)
    if (!charlier::make_path(filtN_path)) {
      charlier::error("filtN_path not created: %s", filtN_path)
      return(RETURN + 1)  
    }
    
    filtN_r <- dadautils::filter_and_trim(fq_files_no_primer, 
                                          output_path  = filtN_path,
                                          maxN         = CFG$dada2_filterAndTrim$maxN, 
                                          multithread  = CFG$multithread, 
                                          truncLen     = CFG$dada2_filterAndTrim$truncLen,
                                          cutoff_params = CFG$dada2_filterAndTrim$cutoff_params,
                                          minLen       = CFG$dada2_filterAndTrim$minLen,
                                          maxLen       = CFG$dada2_filterAndTrim$maxLen,
                                          maxEE        = CFG$dada2_filterAndTrim$maxEE,
                                          truncQ       = CFG$dada2_filterAndTrim$truncQ,
                                          minQ       = CFG$dada2_filterAndTrim$minQ,
                                          rm.phix      = CFG$dada2_filterAndTrim$rm.phix,
                                          compress     = CFG$dada2_filterAndTrim$compress,
                                          verbose      = CFG$dada2_filterAndTrim$verbose,
                                          save_results = FALSE) %>%
      readr::write_csv(file.path(CFG$output_path, "filter_and_trim.csv"))
    
    filtN_files <- auntie::list_filepairs(filtN_path,pattern_forward = "^.*.fastq.gz",verify=F) 
    if(identical(unname(lengths(filtN_files)), c(0,0))){
      stop("No filtN files produced")
    }
    if (!identical(lengths(input_files), lengths(filtN_files))){
      # presumably filter_and_trim dropped some files if we get here... 
      # so we need to trim input_files to match.  We assume that the output basenames are the same as the
      # input basenames - so all we need to do is match
      input_files <- sapply(names(input_files),
                            function(name){
                              ix <- basename(input_files[[name]]) %in% basename(filtN_files[[name]])
                              input_files[[name]][ix]
                            }, simplify = FALSE)
      filtN_r <- filtN_r %>%
        dplyr::filter(reads.out > 0)
      sample.names <- dadautils::extract_sample_names(input_files, rule="basename")
      
      
    } # check for dropped inputs
    
    if ("dada2_learnErrors" %in% names(CFG)){
      charlier::info("learn errors")
      err <- dadautils::learn_errors(filtN_files,
                                     output_path = CFG$output_path,
                                     multithread = CFG$multithread,
                                     save_output = FALSE, 
                                     save_graphics = TRUE,
                                     errorEstimationFunction = PacBioErrfun,
                                     BAND_SIZE = 32) %>%
        dadautils::write_errors(file.path(CFG$output_path, "learn_errors"))
    } # learn errors
  } # filter and trim
  
  CFG$stage <- "supervision"
  CFG$input_path <- filtN_path
  
  charlier::write_config(CFG, file.path(CFG$output_path, basename(CFGFILE)))
  charlier::info("done: %s", CFG$output_path)   
  return(RETURN)
} #main


# we only run is run as a script - not if interactive
if (!interactive()){
  CFGFILE <- commandArgs(trailingOnly = TRUE)[1]
} else {
  CFGFILE <- ""
}

CFG <- charlier::read_config(CFGFILE[1], 
                             autopopulate = TRUE,
                             fields = list(
                               data_path = "data_path",
                               reference_path = "reference_path"),
                             rootname = "global")

if (!interactive()){
  ok <- main(CFG)
  quit(save = "no", status = ok)
}
