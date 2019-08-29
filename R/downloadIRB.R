#' Download IRB datasets
#'
#' Downloads the required IRB datasets.
#' @param output_dir Output directory for the datasets.
#' @param file_dir Path to the compressed IRB database file.
#' @importFrom utils download.file untar
#' @export
downloadIRB <- function(output_dir,
                        file_dir="http://gattaca.imppc.org/genome_browser/lplab/IRB_database.tar.gz") {
  tf <- tempfile()
  message("Will begin downloading datasets to ", tf)
  ret <- download.file(file_dir, tf, mode='wb')

  if (ret != 0) stop("Couldn't download file from ", loc)

  untar(tf, exdir=output_dir, verbose=TRUE)
  message("Done writing IRB dataset files to ", output_dir)

  return(invisible(paste0(output_dir, "IRB_database/")))
}
