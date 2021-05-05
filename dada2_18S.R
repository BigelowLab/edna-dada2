library(dplyr)
library(readr)
library(configr)
library(futile.logger)
library(patchwork)

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
  
  # preliminaries
  charlier::start_logger(filename = file.path(CFG$output_path, "log"))
  PBS_JOBID <- charlier::get_pbs_jobid(no_pbs_text = "not in PBS queue")
  charlier::audit(file.path(CFG$output_path, "audit.txt"), pbs_jobid = PBS_JOBID)
  
  # add baseline info into log, just because
  flog.info("starting run: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  flog.info("NCPUS: %s", as.character(CFG$multithread))
  flog.info("System PID: %s", Sys.getpid())
  flog.info("PBS_JOBID: %s", PBS_JOBID)
  flog.info("VERSION: %s", CFG$version)
  flog.info("INPUT PATH: %s", CFG$input_path)
  flog.info("OUTPUT PATH: %s", CFG$output_path)
  
  flog.info("checking for input fastq files")
  fq_files <- dadautils::list_filepairs(CFG$input_path) %>%
    dadautils::verify_filepairs()
  
  sample.names <- dadautils::extract_sample_names(fq_files, rule="basename")
  
  if ("dada2_plotQualityProfile" %in% names(CFG)){ 
    ofile = file.path(CFG$output_path, "quality_profiles.pdf")
    flog.info("plotting quality profiles: %s", ofile)
    dadautils::plot_qualityProfiles(fq_files, 
                                    n = CFG$dada2_plotQualityProfile$nplots,
                                    ofile = ofile)
  }
  
  flog.info("filter and trim of input files")
  filtN_path <- file.path(CFG$output_path, CFG$dada2_filterAndTrim_filtN$name)
  if (!charlier::make_path(filtN_path)) {
    flog.error("filtN_path not created: %s", filtN_path)
    return(RETURN + 1)  
  }
  
  filtN_r <- dadautils::filter_and_trim(fq_files, 
                                        output_path  = filtN_path,
                                        maxN         = CFG$dada2_filterAndTrim_filtN$maxN, 
                                        multithread  = CFG$multithread, 
                                        truncLen     = CFG$dada2_filterAndTrim_filtN$truncLen,
                                        minLen       = CFG$dada2_filterAndTrim_filtN$minLen,
                                        maxEE        = CFG$dada2_filterAndTrim_filtN$maxEE,
                                        truncQ       = CFG$dada2_filterAndTrim_filtN$truncQ,
                                        rm.phix      = CFG$dada2_filterAndTrim_filtN$rm.phix,
                                        compress     = CFG$dada2_filterAndTrim_filtN$compress,
                                        save_results = TRUE)
  
  
  flog.info("learn errors")
  filtN_files <- dadautils::list_filepairs(filtN_path) 
  
  learnErrors_path <- file.path(CFG$output_path, CFG$dada2_learnErrors$name)
  if (!charlier::make_path(learnErrors_path)) {
    flog.error("learnErrors_path not created: %s", learnErrors_path)
    return(RETURN + 1)  
  }
  
  err <- dadautils::learn_errors(filtN_files,
                                 output_path = learnErrors_path,
                                 multithread = CFG$multithread,
                                 save_output = TRUE, 
                                 save_graphics = TRUE)
  
  
  flog.info("run dada")
  dada_r <- dadautils::run_dada(
    filtN_files, 
    err, 
    multithread = CFG$multithread,
    pool = CFG$dada2_dada$pool)
  
  # run merge pairs
  flog.info("merge pairs")
  mergers <- dadautils::merge_pairs(filtN_files, dada_r, verbose = TRUE, save_output = TRUE, minOverlap = CFG$dada2_merge_pairs$minOverlap)
  #saveRDS(mergers, file = file.path(CFG$output_path, "mergers.rds"))  
  
  
  flog.info("make sequence table")
  seqtab <- dada2::makeSequenceTable(mergers) 
  tseqtab <- dplyr::as_tibble(t(seqtab)) %>%
    readr::write_csv(file.path(CFG$output_path, "seqtab.csv"))
  
  flog.info("remove Bimera Denovo")
  seqtab.nochim <- dada2::removeBimeraDenovo(seqtab, 
                                             method      = CFG$dada2_removeBimeraDenovo_seqtab$method, 
                                             multithread = CFG$multithread, 
                                             verbose     = CFG$dada2_removeBimeraDenovo_seqtab$verbose)
  write.csv(seqtab.nochim, file.path(CFG$output_path, "dada2_seqtab-nochim.csv"))
#  tseqtab.nochim <- dplyr::as_tibble(t(seqtab.nochim)) %>%
#    readr::write_csv(file.path(CFG$output_path, "seqtab-nochim.csv"))
  
  fasta <- dadautils::asv_fasta(seqtab.nochim, file = file.path(CFG$output_path,"ASV_sequences.fasta"))
  
  tseqtab.nochim <- dplyr::as_tibble(t(seqtab.nochim)) %>%
    dplyr::mutate(ASV = names(fasta)) %>%
    dplyr::relocate(ASV, .before = 1) %>%
    readr::write_csv(file.path(CFG$output_path, "seqtab-nochim.csv"))
    
  track <- dplyr::tibble(
    name               = sample.names,
    input              = filtN_r$reads.in, 
    filtered           = filtN_r$reads.out,
    denoised_forward   = sapply(dada_r$forward, dadautils::count_uniques), 
    denoised_reverse   = sapply(dada_r$reverse, dadautils::count_uniques), 
    merged             = sapply(mergers, dadautils::count_uniques), 
    nonchim            = rowSums(seqtab.nochim),
    final_prr          = nonchim/input) %>%
    readr::write_csv(file.path(CFG$output_path, "track.csv"))
  
  
  flog.info("assign taxonomy")
#  taxa <- dadautils::assign_taxonomy(seqtab.nochim, 
#                                     refFasta          = CFG$dada2_assignTaxonomy_nochim$refFasta, 
#                                     taxLevels         = CFG$dada2_assignTaxonomy_nochim$taxLevels, 
#                                     minBoot           = CFG$dada2_assignTaxonomy_nochim$minBoot, 
#                                     outputBootstraps  = CFG$dada2_assignTaxonomy_nochim$outputBootstraps, 
#                                     verbose           = CFG$dada2_assignTaxonomy_nochim$verbose, 
#                                     multithread       = CFG$multithread,
#                                     drop_levels       = "Species")
  
  taxa <- dadautils::assign_taxonomy(seqtab.nochim, 
                                     refFasta          = CFG$dada2_assignTaxonomy_nochim$refFasta, 
                                     taxLevels         = CFG$dada2_assignTaxonomy_nochim$taxLevels, 
                                     minBoot           = CFG$dada2_assignTaxonomy_nochim$minBoot, 
                                     outputBootstraps  = CFG$dada2_assignTaxonomy_nochim$outputBootstraps, 
                                     verbose           = CFG$dada2_assignTaxonomy_nochim$verbose, 
                                     multithread       = CFG$multithread,
                                     drop_levels       = "NA",
                                     save_file         = TRUE,
                                     filename          = file.path(CFG$output_path, "taxa.csv"))

  ttaxa <- dplyr::as_tibble(taxa) %>%
    dplyr::mutate(ASV = names(fasta)) %>%
    dplyr::relocate(ASV, .before = 1) %>%
    readr::write_csv(file.path(CFG$output_path, "ASV_taxa.csv"))

  
  if ("dada2_addSpecies" %in% names(CFG)){
    flog.info("add species to taxonomy")
    if (length(taxa) == 2 && "tax" %in% names(taxa)){
      taxa <- taxa$tax
    }
    taxa <- dada2::addSpecies(taxa, refFasta = CFG$dada2_addSpecies$refFasta)
    readr::write_csv(taxa %>% dplyr::as_tibble(), 
                     file.path(CFG$output_path, "taxa-species.csv"))
  }
  
  if ("dada2_taxa_remove" %in% names(CFG)){
    flog.info("remove unwanted values in taxonomy")
    taxa <- taxa %>%
      dadautils::taxa_remove(vars = CFG$dada2_taxa_remove)
    readr::write_csv(taxa %>% dplyr::as_tibble(), 
                     file.path(CFG$output_path, "taxa-cleaned.csv"))
  }
  
  flog.info("done: %s", basename(filename))  
  return(RETURN)
} #main


# we only run is run as a script - not if interactive
if (!interactive()){
  cfgfile <- commandArgs(trailingOnly = TRUE)[1]
} else {
  cfgfile = ""
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
