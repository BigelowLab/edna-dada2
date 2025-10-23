suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  
  library(ShortRead)  
  library(Biostrings)
  library(dada2)
  
  library(charlier)
  library(dadautils)
})


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
  input_files <- auntie::list_filepairs(CFG$input_path) |>
                   auntie::verify_filepairs()
  if (all( auntie::count_filepairs(input_files) == 0 )) {
    msg <- sprintf("check list_filepairs patterns, no files found in %s", CFG$input_path)
    charlier::error(msg)
    stop(msg)
  }

  # pacbio?
  norev <- auntie::is_singleended(input_file) #(length(input_files) == 1) || (lengths(input_files)[[2]] == 0)
  if (norev) charlier::info("This is a single-ended sample")
  sample.names <- dadautils::extract_sample_names(input_files, rule="before first _")
  
  
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
    if (all(ok == 0)) {
      
      #input_files <- cut_files |>
      #  dadautils::verify_filepairs()
      
      charlier::info("checking for newly created input fastq files within cutadapt path")
      input_files <- auntie::list_filepairs(cutadapt_path) |>
                       auntie::verify_filepairs(min_size=100)
      sample.names <- dadautils::extract_sample_names(input_files,rule="before first _")
      if (all(auntie::count_filepairs(input_files) == 0)) {
        charlier::error("no files made it past cutadapt, check cutadapt settings")
        stop("no files made it past cutadapt, check cutadapt settings")
      }      
    } else {
      charlier::error("Something is wrong with cutadapt: %s", paste(ok, sep = ", "))
      stop("dadutils::run_cutadapt - what just happened to cutadapt?")
    }
  }
  
  
  charlier::info("generate quality profiles")  
  qpp <- quality_profile_pairs(input_files,
                               amplicon_length=CFG$quality$amplicon_length,
                               min_overlap=CFG$quality$min_overlap,
                               n = CFG$quality$sample_n, 
                               params = CFG$dada2_filterAndTrim$cutoff_params,
                               plot_filename=file.path(CFG$output_path, "quality_profiles.pdf"), 
                               overlap_filename=file.path(CFG$output_path, "overlap.csv")) |>
                               dadautils::write_QPP(file.path(CFG$output_path, "qpp"))
  
  # if the user provided explicit truncLen values we use those (this should be very rare or non-existent), if not then we extract those from qpp
  if (inherits(CFG$dada2_filterAndTrim$truncLen, "character")){
    if(CFG$dada2_filterAndTrim$truncLen[1] == "auto"){
      charlier::info("using autocomputed truncLen")
      truncLen <- dplyr::tibble(forward = qpp$forward$cutoff$Cycle,
                                forward_file = input_files$forward)
      if (!norev) truncLen <- dplyr::mutate(truncLen, 
          reverse = qpp$reverse$cutoff$Cycle,
          reverse_file = input_files$reverse)
      trunLen_file <- file.path(CFG$output_path, "truncLen.csv")
      truncLen <- readr::write_csv(truncLen, trunLen_file )
      CFG$dada2_filterAndTrim$truncLen <- trunLen_file 
    } else if (file.exists(CFG$dada2_filterAndTrim$truncLen)){
      charlier::info("reading truncLen from file: %s", CFG$dada2_filterAndTrim$truncLen)
      truncLen <- readr::read_csv(CFG$dada2_filterAndTrim$truncLen, show_col_types = FALSE)
    } else {
      charlier::error("truncLen can be character, bit if so must be 'auto' or a filename")
      stop("truncLen not found:", CFG$dada2_filterAndTrim$truncLen)
    }
  } else {
    charlier::info("using user provided truncLen vector: %s", paste(CFG$dada2_filterAndTrim$truncLen, sep = ", ") )
    truncLen <- CFG$dada2_filterAndTrim$truncLen
  }

  charlier::info("estimate expected error thresholds (for maxEE)")
  peet <- paired_quality_scores(input_files) |>
   paired_ee_per_read() |>
    paired_ee_threshold(sample_names  = sample.names, 
                        filename = file.path(CFG$output_path, "EE_thresholds.csv"))
  
  
  if("dada2_filterAndTrim" %in% names(CFG)){
  
    charlier::info("filter and trim of input files")
    
    filtN_path <- file.path(CFG$output_path, CFG$dada2_filterAndTrim$name)
    if (!charlier::make_path(filtN_path)) {
      charlier::error("filtN_path not created: %s", filtN_path)
      return(RETURN + 1)  
    }
    
    filtN_r <- dadautils::filter_and_trim(input_files, 
                                          output_path  = filtN_path,
                                          maxN         = CFG$dada2_filterAndTrim$maxN, 
                                          multithread  = CFG$multithread, 
                                          truncLen     = truncLen,
                                          cutoff_params = CFG$dada2_filterAndTrim$cutoff_params,
                                          minLen       = CFG$dada2_filterAndTrim$minLen,
                                          maxEE        = CFG$dada2_filterAndTrim$maxEE,
                                          truncQ       = CFG$dada2_filterAndTrim$truncQ,
                                          rm.phix      = CFG$dada2_filterAndTrim$rm.phix,
                                          compress     = CFG$dada2_filterAndTrim$compress,
                                          verbose      = CFG$dada2_filterAndTrim$verbose,
                                          save_results = FALSE) |>
                readr::write_csv(file.path(CFG$output_path, "filter_and_trim.csv"))
    
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
      filtN_r <- filtN_r |>
             dplyr::filter(reads.out > 0)
      sample.names <- dadautils::extract_sample_names(input_files, rule="before first _")
      
      
    } # check for dropped inputs
  
    if ("dada2_learnErrors" %in% names(CFG)){
      charlier::info("learn errors")
      err <- dadautils::learn_errors(filtN_files,
                                     output_path = CFG$output_path,
                                     multithread = CFG$multithread,
                                     save_output = FALSE, 
                                     save_graphics = TRUE) |>
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
  CFGFILE <- "/mnt/storage/data/edna/dada/projects/3step/18S/asv_18S_preprocess.yaml"
}

message(sprintf("CFGFILE: %s", CFGFILE))

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
