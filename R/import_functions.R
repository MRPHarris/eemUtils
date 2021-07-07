# Import functions for use with fluorometer data.

#' An import function for use with processed excitation-emission matrix (PEM) ASCII files from the HORIBA Aqualog. Essentially identical to eemR::eem_import_aqualog, save for a single line.
#'
#' @description To be used within eemR::eem_read() to import aqualog PEM ASCII data files.
#'
#' @param file path to file
#'
#' @export
#'
eem_read_aqualog_PEM <- function(file){

  data <- readLines(file)
  eem <- stringr::str_extract_all(data, "-?\\d+(?:\\.\\d*)?(?:[eE][+\\-]?\\d+)?")
  ex <- sort(as.numeric(eem[[1]]))

  n_col <- lapply(eem, length)
  n_col <- unlist(n_col)
  expected_col <- as.numeric(names(sort(-table(n_col)))[1])

  eem[n_col != expected_col] <- NULL
  eem <- lapply(eem, as.numeric)
  eem <- do.call(rbind, eem)

  em <- eem[, 1]
  eem <- eem[, -1]
  # The line below replaces the hashed line after it.
  eem <- as.matrix(eem) # new as.matrix() without reversal.
  # Original as.matrix() call from eem_import_aqualog(). Led to incorrect reversal of the eem data along the excitation axis.
  #eem <- as.matrix(eem[, ncol(eem):1])

  l <- list(
    file = file,
    x = eem,
    em = em,
    ex = ex
  )
  return(l)
}
