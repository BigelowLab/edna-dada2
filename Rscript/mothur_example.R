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

library(dadautils)


#' main processing step - tries to operate as a pipeline returning 0 (success) or
#' failure ( > 0)
#' 
#' @param cfg character filename of the YAML configuration file
#' @return integer where 0 means success
main <- function(
  cfg = "/home/btupper/edna/edna-dada2/config/mothur_example.yml"){
  
  RETURN <- 0
  CFG <- dadautils::get_configuration(cfg)
  
  if (!dir.exists(CFG$input_path)){
    warning("input path not found:", CFG$input_path)
    return(RETURN + 1)
  }
  
  if (!make_path(CFG$output_path)){
    warning("output path not created:", CFG$output_path)
    return(RETURN + 1)
  }
  
  PBS_JOBID <- Sys.getenv("PBS_JOBID")
  if (nchar(PBS_JOBID) == 0) PBS_JOBID <- "not in PBS queue"
  ok <- dadautils::audit(file.path(CFG$output_path, "audit.txt"), pbs_jobid = PBS_JOBID)
  ok <- flog.threshold(toupper(CFG$verbose[1]))
  ok <- flog.appender(appender.tee(file.path(CFG$output_path, "log")) )
  flog.info("starting run: %s", cfg)
  flog.info("PBS_JOBID: %s", PBS_JOBID)
  flog.info("VERSION: %s", CFG$version)
  flog.info("INPUT PATH: %s", CFG$input_path)
  flog.info("OUTPUT PATH: %s", CFG$output_path)
  
  if (is.numeric(CFG$multithread)){
    MAX_CORES <- dadautils::count_cores()
	  if (interactive() && MAX_CORES < CFG$multithread){
	  	flog.warn("fewer cores available (%i) than requested (%i), adjusting request",
	  			MAX_CORES, CFG$multithread)
	  	CFG$multithread <- MAX_CORES
	  }
  	flog.info("Using %i (of %i) cores", CFG$multithread, MAX_CORES)
  }
  
  if (("cutadapt" %in% names(CFG)) && (nchar(Sys.which(CFG$cutadapt$app)) == 0)){
    flog.error("cutadapt application not found: %s", CFG$cutadapt$app)
    return(RETURN + 1)
  }
  
  flog.info("checking for input fastq files")
  fq_files <- dadautils::list_fastq(CFG$input_path)
  
  if (length(fq_files[[1]]) <= 0){
    flog.error("fastq files not found: %s", CFG$input_path)  
    return(RETURN + 1)
  }
  if (length(fq_files[[1]]) != length(fq_files[[2]])){
    flog.error("unequal number of fastq files: %i forward and %i reverse", 
               length(fq_files$forward), length(fq_files$reverse))
    return(RETURN + 1) 
  }
  sample.names <- sapply(strsplit(basename(fq_files$forward), "_"), `[`, 1)

  ofile <- file.path(CFG$output_path, "quality_profiles.pdf")
  flog.info("plotting quality profiles: %s", ofile)
  ix <- seq_len(CFG$dada2_plotQualityProfile$nplots)
  pdf(ofile)
  	print(dadautils::plot_qualityProfile(fq_files$forward[ix]))
  	print(dadautils::plot_qualityProfile(fq_files$reverse[ix]))
  dev.off()
  

  flog.info("filter and trim of input files")
  filtN_path <- file.path(CFG$output_path, CFG$dada2_filterAndTrim_filtN$name)
  if (!make_path(filtN_path)) {
    flog.error("filtN_path not created: %s", filtN_path)
    return(RETURN + 1)  
  }
  
  filtN_r <- dadautils::filter_and_trim(fq_files, 
                         output_path = filtN_path,
                         maxN        = CFG$dada2_filterAndTrim_filtN$maxN, 
                         multithread = CFG$multithread, 
                         truncLen     = CFG$dada2_filterAndTrim_filtN$truncLen,
                         maxEE       = CFG$dada2_filterAndTrim_filtN$maxEE,
                         truncQ       = CFG$dada2_filterAndTrim_filtN$truncQ,
                         rm.phix      = CFG$dada2_filterAndTrim_filtN$rm.phix,
                         compress    = CFG$dada2_filterAndTrim_filtN$compress) %>%
    readr::write_csv(file.path(CFG$output_path, paste0(CFG$dada2_filterAndTrim_filtN$name, "-results.csv")))
  
  flog.info("learn errors")
  filtN_files <- dadautils::list_fastq(filtN_path)	
  learnErrors_path <- file.path(CFG$output_path$dada2_learnErrors$name)
  if (!make_path(learnErrors_path)) {
    flog.error("learnErrors_path not created: %s", learnErrors_path)
    return(RETURN + 1)  
  }
  err <- dadautils::learn_errors(filtN_files,
  	output_path = learnErrors_path,
  	save_output = TRUE, 
  	save_graphics = TRUE)
  
  flog.info("run dada")
  dada_r <- dadautils::run_dada(
    filtN_files, 
    errs, 
    multithread = CFG$multithread)
  
  # run merge pairs
  flog.info("merge pairs")
  mergers <- dadautils::merge_pairs(filtN_files, dada_r, verbose = TRUE)
  saveRDS(mergers, file = file.path(CFF$output_path, "mergers.rds"))  
  

  flog.info("make sequence table")
  seqtab <- dada2::makeSequenceTable(mergers) 
  tseqtab <- dplyr::as_tibble(t(seqtab)) %>%
    readr::write_csv(file.path(CFG$output_path, "seqtab.csv"))
    
  flog.info("remove Bimera Denovo")
  seqtab.nochim <- dada2::removeBimeraDenovo(seqtab, 
    method      = CFG$dada2_removeBimeraDenovo_seqtab$method, 
    multithread = CFG$multithread, 
    verbose     = CFG$dada2_removeBimeraDenovo_seqtab$verbose)
  tseqtab.nochim <- dplyr::as_tibble(t(seqtab.nochim)) %>%
    readr::write_csv(file.path(CFG$output_path, "seqtab-nochim.csv"))

  track <- dplyr::tibble(
                         name               = sample.names,
                         input              = filtered_r$reads.in, 
                         filtered           = filtered_r$reads.out,
                         denoised_forward   = sapply(dada_r$forward, count_uniques), 
                         denoised_reverse   = sapply(dada_r$reverse, count_uniques), 
                         merged             = sapply(mergers, count_uniques), 
                         nonchim            = rowSums(seqtab.nochim)) %>%
    readr::write_csv(file.path(CFG$output_path, "track.csv"))
  
  flog.info("assign taxonomy")
  taxa <- dada2::assignTaxonomy(seqtab.nochim, 
      refFasta          = CFG$dada2_assignTaxonomy_nochim$refFasta, 
      taxLevels         = CFG$dada2_assignTaxonomy_nochim$taxLevels, 
      minBoot           = CFG$dada2_assignTaxonomy_nochim$minBoot, 
      outputBootstraps  = CFG$dada2_assignTaxonomy_nochim$outputBootstraps, 
      verbose           = CFG$dada2_assignTaxonomy_nochim$verbose, 
      multithread       = CFG$multithread) %>%
    dplyr::as_tibble() %>%
    read::write_csv(file.path(CFG$output_path, "taxa.csv"))

  return(RETURN)
} #main


# we only run is run as a script - not if interactive
if (!interactive()){
  cfgfile <- commandArgs(trailingOnly = TRUE)[1]
  ok <- main(cfgfile)
  quit(save = "no", status = ok)
}
