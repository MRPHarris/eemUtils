#' Get the max intensity value's ex/em position.
#'
#' @description Get the excitation/emission coordinates of the point of maximum fluorescence
#'      intensity within the target EEM.
#'
#' @param eem an eem object compliant with the staRdom/eemR framework
#' @param verbose TRUE/FALSE to return a message with the identified max position.
#'
#' @export
#'
get_eem_max_coords <- function(eem, verbose = FALSE){
  eem_df <- as.data.frame(eem)
  # row of max value
  maxval <- max(eem_df$value, na.rm = TRUE)
  maxrow <- as.numeric(as.character(which.max(eem_df$value)))
  max_ex <- as.numeric(as.character(eem_df$ex[maxrow]))
  max_em <- as.numeric(as.character(eem_df$em[maxrow]))
  maxrowvals <- data.frame(matrix(NA,nrow = 1, ncol = 3))
  colnames(maxrowvals) <- c("em","ex","value")
  maxrowvals[1,1:3] <- c(max_em,max_ex,maxval)
  #maxrowvals <- eem_df[maxrow,1:3]
  if(isTRUE(verbose)){
    message("Max Intensity at ex",max_ex," em",max_em," | value = ",maxval)
  }
  return(maxrowvals)
}
