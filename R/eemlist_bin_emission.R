#' Bin adjacent emission increments/pixels within all EEMs in an eemlist
#'
#' @description Average adjacent emission pixels in EEMs on a per-excitation wavelength basis. Useful for denoising in high-resolution scans.
#'
#' @param eemlist A list of EEMs.
#' @param pix_bins number of adjacent emission increments (pixels on a CCD fluoromter) to group together.
#' @param verbose logical; TRUE/FALSE to print progress. This can be nice as this is not the most efficient set of nested functions - it can take a little while!
#'
#' @export
#'
eemlist_bin_emission <- function(eemlist, pix_bins = 2, verbose = TRUE){
  # Bin em
  eemlist_binned <- vector('list', length = length(eemlist))
  for(e in seq_along(eemlist_binned)){
    if(verbose){message("Binning EEM ",e,"/",length(eemlist))}
    eemlist_binned[[e]] <- eem_bin_emission(eemlist[[e]])
  }
  # eemlist_binned <- lapply(eemlist, function(x){
  #   eem_bin_emission(x, pix_bins = pix_bins)
  # })
  eemlist_binned %>% 'class<-'(c('eemlist'))
  return(eemlist_binned)
}
