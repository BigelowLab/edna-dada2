#' Usage:
#' Rscript ncpus /path/to/output
#'
#' Arguments:
#' ncpus numeric, required, number of CPUs
#' output string, output path, defaults to "." where PBS JOBID stamped subdirectory is made


# capture info about the state
audit <- function(filename = "", pbs_jobid = "not known", ncores = NA){
  cat("[audit]\n")
	cat("Audit date:", 
	    format(Sys.time(), "%Y-%m-%d %H:%M:%S", usetz = TRUE), "\n",
		   file = filename, append = TRUE)
	cat("System PID:", Sys.getpid(), "\n", file = filename, append = TRUE)
	cat("PBS_JOBID:", pbs_jobid, "\n", file = filename, append = TRUE)
	cat("Cores:", ncores, "\n", file = filename, append = TRUE)
	cat("R version:", R.version.string, "\n", file = filename, append = TRUE)
	cat("libPaths():\n", file = filename, append = TRUE)
	for (lp in .libPaths()) {
	  cat("    ", lp, "\n", file = filename, append = TRUE)
	}
	x <- as.data.frame(utils::installed.packages(), stringsAsFactors = FALSE)
	x <- x[,c("Package", "Version",  "LibPath")]
	cat("[installed.packages]\n", file = filename, append = TRUE)
	if (nzchar(filename)){
		conn <- file(filename, open = 'at')
		utils::write.csv(x, file = conn, row.names = FALSE)
		close(conn)
	} else {
		print(x, row.names = FALSE)
	}
	cat("[end-audit]\n", file = filename, append = TRUE)
	invisible(NULL)
}

#simplified message output
catn <- function(..., file = LOGFILE, append = TRUE){
  cat(..., "\n", file = file, append = append)
}


library(dada2)

# write all output to file
JOBID <- Sys.getenv("PBS_JOBID")
if (nchar(JOBID) <= 0) JOBID <- format(Sys.time(), "%Y%m%d%H%m%s") 

# read arguments from command line
args = commandArgs(trailingOnly=TRUE)
NCORES <- as.numeric(args[1])
if (length(args) > 1){
  OPATH <- args[2]
} else {
  OPATH <- "."
} 
LOGFILE <- file.path(OPATH, paste0(JOBID, ".Rout"))
audit(filename = LOGFILE, pbs_jobid = JOBID, ncores = NCORES)

IPATH <- "/mnt/storage/labs/countway_nfs/dada2_data/SK18S_batch/cutadapt/"
catn(sprintf("IPATH = %s", IPATH))

catn(sprintf("listing original input files in %s", IPATH))
fnFs <- sort(list.files(IPATH, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(IPATH, pattern="_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
catn(sprintf("original sample.names = %s", paste(sample.names, collapse = " ")))

ipath <- file.path(IPATH, "filtered")
catn(sprintf("listing original filtered files in %s", ipath))
filtFs <- file.path(ipath, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(ipath, paste0(sample.names, "_R_filt.fastq.gz"))
filtF.sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
filtR.sample.names <- sapply(strsplit(basename(filtRs), "_"), `[`, 1)
catn(sprintf("filtered forward sample.names = %s", 
             paste(filtF.sample.names, collapse = " ")))
catn(sprintf("filtered reverse sample.names = %s", 
             paste(filtR.sample.names, collapse = " ")))


catn("learn errors forward...")
sink(LOGFILE, append = TRUE)
errF <- learnErrors(filtFs, multithread=NCORES)
sink(NULL)

catn("learn errors reverse...")
sink(LOGFILE, append = TRUE)
errR <- learnErrors(filtRs, multithread=NCORES)
sink(NULL)

catn("learn errors complete...")