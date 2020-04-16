#' Count the number of matches of the primer in the subject
#' 
#' @export
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
#' @export
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
#' @export
#' @param primer input sequence, see \code{\link[Biostrings]{DNAString}}
#' @return name character vector
all_orients <- function(primer = "TTGAAAA-CTC-N")  {
  dna <- Biostrings::DNAString(primer)
  orients <- c(Forward = dna, 
             Complement =  Biostrings::complement(dna), 
             Reverse = IRanges::reverse(dna), 
             RevComp = Biostrings::reverseComplement(dna)) 
  return(sapply(orients, toString))
} 