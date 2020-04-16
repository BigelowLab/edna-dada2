#' Run cutadapt
#'
#' @export
#' @param cut_files two element list for forward and reverse cutadapt results to be created
#' @param filt_files two element list of forward and reverse filtered files
#' @param CFG list of configuration
#' @param save_output logical, if TRUE try to capture the output of cutadapt to text files,
#'   one per cutadapt file generated, in the cutadapt output directory
#' @param save_graphics logical, if TRUE try to capture quality plots form the resulting cut_files
#' @return numeric codes, one per cut_file pairing as per \code{system2} where 0 means success
run_cutadapt <- function(
  cut_files,
  filt_files,
  CFG,
  save_output = TRUE,
  save_graphics = FALSE){
  
  FWD.RC <- dada2::rc(CFG$primer$FWD)
  REV.RC <- dada2::rc(CFG$primer$REV)
  R1.flags <- paste("-g", CFG$primer$FWD, "-a", REV.RC)
  R2.flags <- paste("-G", CFG$primer$REV, "-A", FWD.RC)

  OK <- sapply(seq_along(cut_files$forward),
    function(i){
      if (save_output){
        ofile <- paste0(strip_extension(cut_files$forward[1]), ".cutadapt_output.txt")
      } else {
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
    ofile <- paste0(strip_extension(cut_files$forward[1]), ".cutadapt_quality.pdf")
    grDevices::pdf(ofile)
    try(dada2::plotQualityProfile(cut_files$forward[ix]) +  dada2::plotQualityProfile(cut_files$reverse[ix]))
    grDevices::dev.off()
  }
  OK
}