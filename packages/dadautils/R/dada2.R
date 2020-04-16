#' Filter and trim
#' 
#' @export
#' @param filelist list of forward and reverse fastq files
#' @param output_path character, the output path
#' @param ... other arguments for \code{\link[dada2]{filterAndTrim}}
#' @param save_results logical, save CSV if TRUE
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


#' Run dada2::learnErrors on a set of fastq files
#'
#' @export
#' @param filelist list of forward and reverse fastq files
#' @param ... arguments for \code{\link[dada2]{learnErrors}}
#' @param output_path character, the output path
#' @param save_output logical, if TRUE save the output to the specified output_path
#' @param save_graphics logical, if TRUE try to capture quality plots form the resulting cut_files
learn_errors <- function(filelist, ...,
  output_path = dirname(filelist$forward[1]),
  save_output = FALSE, 
  save_graphics = FALSE){

  errs <- list(
      forward =  dada2::learnErrors(filelist$forward, ...),
      reverse =  dada2::learnErrors(filelist$reverse, ...)
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
    grDevices::pdf(ofile, height = 7.5, width = 10.5)    
    try(
      print(pforward + preverse)
    )
    grDevices::dev.off()
  }

  errs
}


#' Run dada
#'
#' @export
#' @param filelist list of forward and reverse fastq files
#' @param errs list of forward and reverse outputs of learnErrors
#' @param ... arguments for \code{\link[dada2]{dada}}
#' @return list with elements for forward and reverse as returned by \code{\link[dada2]{dada}}
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
#' @export
#' @param filelist list of forward and reverse fastq files
#' @param dada_r  list of dada2::dada results
#' @param ... arguments for \code{\link[dada2]{mergePairs}}
#' @return as returned by \code{\link[dada2]{mergePairs}}
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


#' Count uniques
#'
#' @export
#' @param x object from which uniques-vector can be extracted
#' @param ... further arguments for  \code{\link[dada2]{getUniques}}
#' @return integer number of uniques
count_uniques <- function(x, ...){
  sum(dada2::getUniques(x, ...)) 
}



#' Plot quality profiles for one or more FASTQs
#'
#' @export
#' @param x character vector of fastq filenames 
#' @param ofile character, the name of the PDF file to generate of NA
#' @param ... further arguments for \code{\link[dada2]{plotQualityProfile}}
plot_qualityProfile <- function(x, 
	ofile = c(NA, "qualityProfile.pdf")[1],
	...){
	if (!is.na(ofile)) grDevices::pdf(ofile)
	ok <- dada2::plotQualityProfile(x, ...)
	if (!is.na(ofile)) grDevices::dev.off()
	ok
}

