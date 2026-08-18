#' A single-point scatter denoiser for EEM measurements.
#'
#' @description EEMs are often punctuated by isolated, single points of scatter. These are
#'      potentially caused by a wide range of factors, including gamma-ray spikes. They can
#'      severely hamper normalisation attempts, and throw off method detection limit calculations.
#'      This function provides a simple denoising solution using the tsclean function from the
#'      forecast package. Data supplied to this function should not contain missing values.
#'
#' @param eemlist A list of EEMs in a format compliant with eemR/staRdom.
#' @param interp TRUE/FALSE to replace outliers with linear interpolated values instead of NA.
#'
#' @export
#'
eemlist_sp_denoise <- function(eemlist,
                               interp = TRUE){
  denoise_eemlist <- lapply(eemlist, function(e){
    e_df <- as.data.frame(e, gather = FALSE)
    clean_df <- lapply(e_df, function(ex_wl){
      wl_clean <- forecast::tsclean(ex_wl, replace.missing = interp)
    }) %>%
      data.frame() %>%
      'rownames<-'(rownames(e_df)) %>%
      'colnames<-'(colnames(e_df))
    clean_df <- eemdf_to_eem(eemdf = clean_df,
                             file = e$file,
                             sample = e$sample,
                             location = e$location,
                             gathered = FALSE)
    clean_df
  })
  return(denoise_eemlist)
}
