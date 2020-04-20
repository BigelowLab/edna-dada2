#' Provide an R session audit
#'
#' @param filename the name of the file to dump to or "" to print to console
audit <- function(filename = "", pbs_jobid = "not known"){
	cat("Audit date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S", usetz = TRUE), "\n",
		file = filename)
	cat("PBS_JOBID:", pbs_jobid, "\n", file = filename, append = TRUE)
	cat("R instance:", Sys.which("R"), "\n", file = filename, append = TRUE)
	cat("R version:", R.version.string, "\n", file = filename, append = TRUE)
	cat("libPaths():\n", file = filename, append = TRUE)
	for (lp in .libPaths()) cat("    ", lp, "\n",
															file = filename, append = TRUE)
	x <- as.data.frame(installed.packages(), stringsAsFactors = FALSE)
	x <- x[,c("Package", "Version",  "LibPath")]
	cat("installed.packages():\n",
			file = filename, append = TRUE)
	if (nzchar(filename)){
		conn <- file(filename, open = 'at')
		write.csv(x, file = conn, row.names = FALSE)
		close(conn)
	} else {
		print(x, row.names = FALSE)
	}
	invisible(NULL)
}


if (!interactive()){

	path <- commandArgs(trailingOnly = TRUE)[1]
	if (length(path) == 0 || nchar(path[1]) == 0 ) {
		path <- "/home/btupper/edna/data/audits"
	}
	PBS_JOBID <- Sys.getenv("PBS_JOBID")
  if (nchar(PBS_JOBID) == 0){
  	PBS_JOBID <- "not in PBS queue"
  	filename <- format(Sys.time(), "%Y-%m-%dT%H%m%s.txt")
  } else {
  	filename <- paste0(PBS_JOBID, ".txt")
  }
	
	audit(filename = file.path(path, filename),
				pbs_jobid = PBS_JOBID)
}