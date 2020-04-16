



#' Count the number of CPUs
#'
#' THis is a wrapper around \code{\link[parallel]{detectCores}}
#'
#' @return integer
count_cores <- function(){
  parallel::detectCores()
}


#' List fastq files and separate into forward and reverse reads
#' 
#' @param path character, the input path
#' @param pattern_forward file pattern 
#' @param pattern_reverse file pattern
#' @return named list of sorted foreward and reverse fastq filenames
list_fastq <- function(path,
                       pattern_forward = "^.*_F_filt\\.fastq",
                       pattern_reverse = "^.*_R_filt\\.fastq"){
  
  list(
    forward = sort(list.files(path, pattern = pattern_forward, full.names = TRUE)),
    reverse = sort(list.files(path, pattern = pattern_reverse, full.names = TRUE)) )
}


do_order <- tolower(commandArgs(trailingOnly = TRUE))
if (length(do_order) == 0){
	do_order = c("forward", "reverse")
}
if (length(do_order) == 1){
	two <- switch(do_order,
		"forward" = "reverse",
		"reverse" = "forward")
	do_order <- c(do_order, two)
}


path <- "/home/btupper/edna/data/examples/SK18S/filtered"
sink( file.path(path, sprintf("learn_errors_%s_%s.out", do_order[1], do_order[2]) ) )

filelist <- list_fastq(path)
verbose <- 2
multithread = count_cores()

cat("do_order:", do_order[1], do_order[2], "\n")
cat("path:", path, "\n")
cat("verbose:", verbose, "\n")
cat("multithread:", multithread, "\n")

for (d in do_order){
	cat("\n\n***", d, "***\n")
	cat(paste("  ", filelist[[d]]), sep = "\n")
	result <- dada2::learnErrors(filelist[[d]], multithread = multithread, verbose = verbose)
}

sink(NULL)
    