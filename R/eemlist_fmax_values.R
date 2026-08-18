#' Obtain the fmax value for each EEM in an eemlist.
#'
#' @description obtain the maximum flurescence intensity value (fmax) from a set of EEMs.
#'
#' @param eemlist a list of eem objects compliant with the staRdom/eemR framework.
#'
#' @export
#'
eemlist_fmax_values <- function(eemlist){
  valtab <- data.frame(matrix(NA,nrow = length(eemlist), ncol = 1))
  for (i in seq_along(eemlist)) {
    sample_name <- eemlist[[i]][["sample"]]
    filename <- eemlist[[i]][["file"]]
    location <- eemlist[[i]][["location"]]
    eem_ungathered <- as.data.frame(eemlist[[i]], gather = FALSE)
    valtab[i,1] <- max(eem_ungathered, na.rm = TRUE)
  }
  valtab
}
