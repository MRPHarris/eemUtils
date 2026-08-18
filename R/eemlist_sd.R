#' Generate standard deviations of an eemlist
#'
#' @description Produces an EEM object containing ex/em-pair standard deviations for all eems within the supplied eemlist.
#'
#' @param eemlist A list of EEMs in a format compliant with eemR/staRdom.
#' @param mult Numeric; a simple numeric multiplier applied to the resulting ex-em pair sd values.
#'
#' @export
#'
eemlist_sd <- function(eemlist, mult = 1){
  cols <- eemlist[[1]]$ex
  rows <- eemlist[[1]]$em
  # Coerce eemlist to sd-operable matrices
  eemlist_matrices <- lapply(eemlist, function(x){
    eem_df <- as.data.frame(x, gather = FALSE)
    eem_df <- as.matrix(eem_df)
  })
  # Produce standard deviations
  sd <- apply(simplify2array(eemlist_matrices), 1:2, sd)
  # Multiplication handling
  if(mult != 1){
    if(!is.numeric(mult)){
      stop("'mult' must be numeric.")
    } else {
      sd <- sd * mult
    }
  }
  sd <- data.frame(sd)
  colnames(sd) <- cols
  rownames(sd) <- rows
  sd
}
