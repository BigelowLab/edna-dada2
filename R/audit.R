#' Provide an R session audit
#'
#' @param filename the name of the file to dump to or "" to print to console
audit <- function(filename = ""){
	cat("Audit date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S", usetz = TRUE), "\n",
		file = filename)
	cat("R version:", R.version.string, "\n",
		file = filename, append = TRUE)
	cat("libPaths():\n",
		file = filename, append = TRUE)
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

