supressPackageStarupMessages({
  library(dplyr)
  library(readr)

  library(ShortRead)  
  library(Biostrings)
  library(dada2)

  library(charlier)
  library(dadautils)
  )}

#' main processing step - tries to operate as a pipeline returning 0 (success) or
#' failure ( > 0)
#' 
#' @param CFG list, configuration from a config file
#' @return integer where 0 means success
main <- function(CFG){
  
  RETURN <- 0
  
  if (tolower(CFG$stage) != "postprocess"){
    warning("configuration stage should be 'postprocess'")
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
  charlier::audit(file.path(CFG$output_path, "audit_postsupervision.txt"), pbs_jobid = PBS_JOBID)
  
  # add baseline info into log, just because
  charlier::info("**** starting run ****")
  charlier::info("NCPUS: %s", as.character(CFG$multithread))
  charlier::info("System PID: %s", Sys.getpid())
  charlier::info("PBS_JOBID: %s", PBS_JOBID)
  charlier::info("VERSION: %s", CFG$version)
  charlier::info("INPUT PATH: %s", CFG$input_path)
  charlier::info("OUTPUT PATH: %s", CFG$output_path)
  
  charlier::info("checking for input fastq files")
  input_files <- dadautils::list_filepairs(CFG$input_path) %>%
                   dadautils::verify_filepairs()

  # pacbio?
  norev <- (length(input_files) == 1) || (lengths(input_files)[[2]] == 0)
 
  sample.names <- dadautils::extract_sample_names(input_files, rule="before first _")
  
  
  filter_trim_file <- file.path(CFG$output_path, "filter_and_trim.csv")
  if (file.exists(filter_trim_file)){
    charlier::info("read filter_and_trim")
    filter_trim <- suppressMessages(readr::read_csv(filter_trim_file))
  } else {
    charlier::error("filter_and_trim file not found: %s", filter_trim_file)
    stop("filter_and_trim file not found: ", filter_trim_file)
  }
  
  
  err_file <- file.path(CFG$output_path, "learn_errors")
  if (file.exists(err_file)){
    charlier::info("read learn errors")
    err <- dadautils::read_errors(err_file)
  } else {
    charlier::error("learn_errors file not found: %s", err_file)
    stop("learn_errors file not found: ", err_file)
  }
  
  charlier::info("run dada")
  dada_r <- dadautils::run_dada(input_files, err, 
                                multithread = CFG$multithread,
                                pool = CFG$dada2_dada$pool)
  
  # run merge pairs
  charlier::info("merge pairs")
  mergers <- dadautils::merge_pairs(input_files, dada_r, 
                                    verbose = TRUE,
                                    minOverlap = CFG$dada2_merge_pairs$minOverlap) %>%
             dadautils::write_mergers(file.path(CFG$output_path, "mergers"))
  
  
  charlier::info("make sequence table")
  seqtab <- dada2::makeSequenceTable(mergers) 
  tseqtab <- dplyr::as_tibble(t(seqtab)) %>%
    readr::write_csv(file.path(CFG$output_path, "seqtab.csv"))
  
  charlier::info("remove Bimera Denovo")
  seqtab.nochim <- dada2::removeBimeraDenovo(seqtab, 
                                             method      = CFG$dada2_removeBimeraDenovo_seqtab$method, 
                                             multithread = CFG$multithread, 
                                             verbose     = CFG$dada2_removeBimeraDenovo_seqtab$verbose)

  
  charlier::info("write fasta")
  fasta <- dadautils::asv_fasta(seqtab.nochim, file = file.path(CFG$output_path,"ASV_sequences.fasta"))
  
  tseqtab.nochim <- dplyr::as_tibble(t(seqtab.nochim)) %>%
    dplyr::mutate(ASV = names(fasta)) %>%
    dplyr::relocate(ASV, .before = 1) %>%
    readr::write_csv(file.path(CFG$output_path, "seqtab-nochim.csv"))
    
  charlier::info("write track")
  track <- dplyr::tibble(
                         name               = sample.names,
                         input              = filter_trim$reads.in, 
                         filtered           = filter_trim$reads.out,
                         denoised_forward   = sapply(dada_r$forward, dadautils::count_uniques), 
                         denoised_reverse   = sapply(dada_r$reverse, dadautils::count_uniques), 
                         merged             = sapply(mergers, dadautils::count_uniques), 
                         nonchim            = rowSums(seqtab.nochim),
                         final_prr          = nonchim/input) %>%
    readr::write_csv(file.path(CFG$output_path, "track.csv"))
  
  
  charlier::info("assign taxonomy")
  # always returns a 2-element list with tax and boot as matrices (or boot = NULL)
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
  
  charlier::info("writing ASV_taxa")
  ttaxa <- dplyr::as_tibble(taxa$tax) %>%
    dplyr::mutate(ASV = names(fasta)) %>%
    dplyr::relocate(ASV, .before = 1) %>%
    dplyr::bind_cols(dplyr::as_tibble(taxa$boot)) %>%
    readr::write_csv(file.path(CFG$output_path, "ASV_taxa.csv"))

  
  if ("dada2_addSpecies" %in% names(CFG)){
    charlier::info("add species to taxonomy")
    if (length(taxa) == 2 && "tax" %in% names(taxa)){
      taxa <- taxa$tax
    }
    taxa <- dada2::addSpecies(taxa, refFasta = CFG$dada2_addSpecies$refFasta)
    readr::write_csv(taxa %>% dplyr::as_tibble(), 
                     file.path(CFG$output_path, "taxa-species.csv"))
  }
  
  if ("dada2_taxa_remove" %in% names(CFG)){
    charlier::info("remove unwanted values in taxonomy")
    taxa <- taxa %>%
      dadautils::taxa_remove(vars = CFG$dada2_taxa_remove)
    readr::write_csv(taxa %>% dplyr::as_tibble(), 
                     file.path(CFG$output_path, "taxa-cleaned.csv"))
  }
  
  
  CFG$stage <- "dada-complete"
  charlier::write_config(CFG, file.path(CFG$output_path, paste0("complete-", basename(CFGFILE))))
  charlier::info("done: %s", CFG$output_path)   
  return(RETURN)
} #main


# we only run is run as a script - not if interactive
if (!interactive()){
  CFGFILE <- commandArgs(trailingOnly = TRUE)[1]
} else {
  CFGFILE <- "/mnt/storage/data/edna/dada/projects/3step/18S/asv_18S_postprocess.yaml"
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
