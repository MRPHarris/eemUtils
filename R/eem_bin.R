#' Bin the intensity values from an EEM.
#'
#' @description Takes the intensity values of an EEM and bins them based upon a preset number of
#'        bins. The lowest break-end is the smallest intensity value, and the highest will be Inf, such
#'        that the top bin is
#'        (max intensity - ((max intensity)/nbins)):Inf
#'
#' @param eem An eem object, compliant with the eemR/staRdom packages.
#' @param nbins The number of bins. 12 by default.
#'
#' @export
#'
eem_bin <- function(eem,
                    nbins = 12){
  if(is(eem,"data.frame")){
    eemdf <- eem
  } else if(is(eem, "eem")){
    eemdf <- as.data.frame(eem, gather = TRUE)
  }
  if(colnames(eemdf)[3] != "value"){
    stop("Please use a gathered EEM dataframe returned by as.data.frame(eem, gather = TRUE)")
  }
  eemdf$value <- as.numeric(eemdf$value)
  max_eem_val <- max(eemdf$value, na.rm = TRUE)
  bin1 <- max_eem_val/(nbins)
  intensity_breaks <- seq(0,max_eem_val-bin1,bin1)
  intensity_breaks <- append(intensity_breaks, Inf, after = length(intensity_breaks))
  if(sum(eem$value < 0, na.rm = TRUE) != 0){
    # handling for negative values
    intensity_breaks <- append(intensity_breaks, min(eemdf$value, na.rm = TRUE), after = 0)
  }
  intensity_labels <- head(intensity_breaks, -1)
  #intensity_labels <- intensity_breaks
  data.table::setDT(eemdf)
  data.table::setDT(eemdf)[ , intensitygroups := cut(value,
                                                     breaks = intensity_breaks,
                                                     right = FALSE,
                                                     labels = intensity_labels)]
  eemdf$value <- as.numeric(as.character(eemdf$intensitygroups))
  eemdf_2 <- subset(eemdf, select = -c(intensitygroups))
  return(eemdf_2)
}
