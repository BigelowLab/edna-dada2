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
  charlier:info("starting run: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  charlier:info("NCPUS: %s", as.character(CFG$multithread))
  charlier:info("System PID: %s", Sys.getpid())
  charlier:info("PBS_JOBID: %s", PBS_JOBID)
  charlier:info("VERSION: %s", CFG$version)
  charlier:info("INPUT PATH: %s", CFG$input_path)
  charlier:info("OUTPUT PATH: %s", CFG$output_path)
  
  charlier:info("checking for input fastq files")
  input_files <- dadautils::list_filepairs(CFG$input_path) %>%
                   dadautils::verify_filepairs()
  
  sample.names <- dadautils::extract_sample_names(input_files, rule="basename")
  
  # cutadapt?
  if ("cutadapt" %in% names(CFG)){
    charlier::info("run cutadapt on input files")
    cutadapt_path <- file.path(CFG$output_path, "cutadapt")
    ok <- charlier::make_path(cutadapt_path)
    cut_files <- sapply(input_files,
      function(ff){
        file.path(cutadapt_path, basename(ff))  
      }, simplify = FALSE)
    ok <- dadautils::run_cutadapt(cut_files, input_files, CFG, save_graphics = FALSE)
    if (ok == 0) {
      input_files <- cut_files
    } else {
      charlie::error("Something is wrong with cutadapt")
      stop("what just happened to cutadapt?")
    }
  }
  
  
  charlier:info("generate quality profiles")  
  qpp <- quality_profile_pairs(input_files, 
                               plot_filename=file.path(CFG$output_path, "quality_profiles.pdf"), 
                               overlap_filename=file.path(CFG$output_path, "overlap.csv"))

  charlier:info("estimate expected error thresholds (for maxEE)")
  pqs <- paired_quality_scores(input_files) %>%
   paired_ee_per_read() %>%
    paired_ee_threshold(sample_names  = sample.names, 
                        filename = file.path(CFG$output_path, "EE_thresholds.csv"))
  
  
  if(CFG$dada2_filterAndTrim_filtN$run){
  
    charlier:info("filter and trim of input files")
    filtN_path <- file.path(CFG$output_path, CFG$dada2_filterAndTrim_filtN$name)
    if (!charlier::make_path(filtN_path)) {
      charlier:error("filtN_path not created: %s", filtN_path)
      return(RETURN + 1)  
    }
    
    filtN_r <- dadautils::filter_and_trim(input_files, 
                                          output_path  = filtN_path,
                                          maxN         = CFG$dada2_filterAndTrim_filtN$maxN, 
                                          multithread  = CFG$multithread, 
                                          truncLen     = CFG$dada2_filterAndTrim_filtN$truncLen,
                                          cutoff_params = CFG$dada2_filterAndTrim_filtN$cutoff_params,
                                          minLen       = CFG$dada2_filterAndTrim_filtN$minLen,
                                          maxEE        = CFG$dada2_filterAndTrim_filtN$maxEE,
                                          truncQ       = CFG$dada2_filterAndTrim_filtN$truncQ,
                                          rm.phix      = CFG$dada2_filterAndTrim_filtN$rm.phix,
                                          compress     = CFG$dada2_filterAndTrim_filtN$compress,
                                          verbose      = CFG$dada2_filterAndTrim_filtN$verbose,
                                          save_results = TRUE)
    
    
    
    filtN_files <- dadautils::list_filepairs(filtN_path) 
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
  
    if (CFG$dada2_learnErrors$run){
      charlier:info("learn errors")
      
      learnErrors_path <- file.path(CFG$output_path, CFG$dada2_learnErrors$name)
      if (!charlier::make_path(learnErrors_path)) {
        charlier:error("learnErrors_path not created: %s", learnErrors_path)
        return(RETURN + 1)  
      }
      
      err <- dadautils::learn_errors(filtN_files,
                                     output_path = learnErrors_path,
                                     multithread = CFG$multithread,
                                     save_output = TRUE, 
                                     save_graphics = TRUE)
    } # learn errors
  } # filter and trim
  
  charlier:info("done: %s", CFG$output_path)   
  return(RETURN)
} #main


# we only run is run as a script - not if interactive
if (!interactive()){
  cfgfile <- commandArgs(trailingOnly = TRUE)[1]
} else {
  cfgfile = "/mnt/storage/data/edna/packages/edna-dada2/avs_18S.yaml"
}

CFG <- charlier::read_config(cfgfile[1], 
                             autopopulate = TRUE,
                             fields = list(
                             data_path = "data_path",
                             reference_path = "reference_path"),
                             rootname = "global")



if (!interactive()){
  ok <- main(CFG)
  quit(save = "no", status = ok)
}
