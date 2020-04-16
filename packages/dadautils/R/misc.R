#' Provide an R session audit
#'
#' @export
#' @param filename character, the name of the file to dump to or "" to print to console
#' @param pbs_jobid character, the OPBS jobid if known
#' @return NULL invisibly
audit <- function(filename = "", pbs_jobid = "not known"){
	cat("Audit date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S", usetz = TRUE), "\n",
		file = filename)
	cat("PBS_JOBID:", pbs_jobid, "\n", file = filename, append = TRUE)
	cat("Cores:", count_cores(), "\n", file = filename, append = TRUE)
	cat("R version:", R.version.string, "\n",
		file = filename, append = TRUE)
	cat("libPaths():\n",
		file = filename, append = TRUE)
	for (lp in .libPaths()) cat("    ", lp, "\n",
															file = filename, append = TRUE)
	x <- as.data.frame(utils::installed.packages(), stringsAsFactors = FALSE)
	x <- x[,c("Package", "Version",  "LibPath")]
	cat("installed.packages():\n",
			file = filename, append = TRUE)
	if (nzchar(filename)){
		conn <- file(filename, open = 'at')
		utils::write.csv(x, file = conn, row.names = FALSE)
		close(conn)
	} else {
		print(x, row.names = FALSE)
	}
	invisible(NULL)
}

#' Count the number of CPUs
#'
#' THis is a wrapper around \code{\link[parallel]{detectCores}}
#'
#' @export
#' @return integer count of cores
count_cores <- function(){
  parallel::detectCores()
}

#' Given a path - make it if it doesn't exist
#' 
#' We don't use R's \code{dir.create} in case we a reaching across volumes which 
#' has caused headaches at time.
#' 
#' @export
#' @param path character, the path to check and/or create
#' @return logical, TRUE if the path exists
make_path <- function(path){
  ok <- dir.exists(path[1])
  if (!ok){
    ok <- system2("mkdir", args = path[1]) == 0
  }
  ok
}


#' Retrieve a default configuration it's just an example
#'
#' @export
#' @return configuration list
default_configuration <- function(){

list(version = "v0.000", 
     email = "btupper@bigelow.org", 
     input_path = "/home/btupper/edna/data/examples/ben_demo_raw", 
     output_path = "/home/btupper/edna/data/examples/ben_demo_raw_results", 
     verbose = "info", 
     multithread = 24L, 
     dada2_plotQualityProfiles = list(
        `FALSE` = 2L), 
     dada2_filterAndTrim_filtN = list(
        name = "filtN", 
        maxN = 0, 
        compress = FALSE), 
    dada2_learnErrors = list(
        name = "learnErrors"), 
    cutadapt = list(
        app = "/mnt/modules/bin/dada2/1.12/bin/cutadapt", 
        more_args = "--minimum-length 25.0 -n 2.0"), 
    dada2_filterAndTrim_filtered = list(
        name = "filtered", 
        truncLen = c(275L, 225L), 
        maxN = 0, 
        maxEE = c(2L, 2L), 
        truncQ = 2L, 
        rm.phix = TRUE, 
        compress = FALSE), 
    dada2_dada_filtered = list(
    	name = "dada"), 
    dada2_removeBimeraDenovo_seqtab = list(
        method = "consensus", 
        verbose = TRUE), 
    dada2_assignTaxonomy_nochim = list(
        refFasta = "/home/btupper/edna/data/examples/pr2_version_4.12.0_18S_dada2.fasta", 
        minBoot = 0L, 
        outputBootstraps = TRUE, 
        verbose = TRUE, 
        taxLevels = c("Kingdom", "Supergroup", "Division", "Class", 
        "Order", "Family", "Genus", "Species")), 
    primer = list(
        FWD = "CYGCGGTAATTCCAGCTC", 
        REV = "AYGGTATCTRATCRTCTTYG"))

}


#' Check the configuration to make sure it is complete, borrowing from defaults as needed.
#'
#' @export
#' @param x list, configuration list
#' @param default, list default configuration
#' @return updated configuration
check_configuration <- function(
  x = list(foo = list(bar = 7, biz = "h"), 
           dog = list(drool = TRUE)), 
  default = list(foo = list(bar = 9, biz = "p"), 
                 fish = list(pickled = "yum", breakfast = TRUE),
                 dog = list(drool = TRUE, quantity = "a lot"),
                 cat = list(hairball = TRUE))){

  nm_d <- names(default)
  ix <- nm_d %in% names(x)
  if (any(!ix))  x[nm_d[!ix]] <- default[nm_d[!ix]]
  
  for (nm in nm_d){
    nms <- names(default[[nm]])
    ix <- nms %in% names(x[[nm]])
    if (any(!ix)){
      wix <- which(ix)
      for (i in wix) x[[nm]][[nms[i]]] <- default[[nm]][[nms[i]]]
    }
  }
  x
}


#' Get the configuration use default if needed
#' 
#' @export
#' @param x configuration filename
#' @param default list describing the default
#' @param check logical, if TRUE check the config against the default
#' @return list config values or a try-error
get_configuration <- function(x ,
  default = default_configuration(),
  check = FALSE){
  
  if (length(x) == 0) {
    cfg <- default
  } else {
    cfg <- try(configr::read.config(x[1]))
    if (inherits(cfg, 'try-error')){
      print(cfg)
      stop("failed to read config file:", x[1])
    }
    # TODO
    if (check) cfg <- check_configuration(cfg, default)
  }
  cfg
}


#' Strip the extension(s) off of a filename
#'
#' Note if ext is ".fastq" then ".fastq.gz" and ".fastq.tar.gz" will also be stripped
#'
#' @export
#' @param filename character one or more filenames
#' @param ext character, one or more extension patterns
#' @return filename with extension stripped
strip_extension <- function(
  filename = c("BR2_2016_S216_L001_R2_001.fastq", "foobar.fastq.gz", "fuzzbaz.txt"),
  ext = ".fastq"){
  
  ix <- gregexpr(ext, filename, fixed = TRUE)
  sapply(seq_along(ix), 
    function(i){
      if (ix[[i]] != -1) {
        s <- substring(filename[i], 1, ix[[i]]-1)
      } else {
        s <- filename[i]
      }
      s
    })
  
  
  }

#' List fastq files and separate into forward and reverse reads
#' 
#' @export
#' @param path character, the input path
#' @param pattern_forward file pattern 
#' @param pattern_reverse file pattern
#' @return named list of sorted foreward and reverse fastq filenames
list_fastq <- function(path,
                       pattern_forward = "^.*_R1_001\\.fastq",
                       pattern_reverse = "^.*_R2_001\\.fastq"){
  
  list(
    forward = sort(list.files(path, pattern = pattern_forward, full.names = TRUE)),
    reverse = sort(list.files(path, pattern = pattern_reverse, full.names = TRUE)) )
}


