#' Bin adjacent emission increments/pixels in an EEM
#'
#' @description Average adjacent emission pixels in an EEM on a per-excitation wavelength basis. Useful for denoising in high-resolution scans.
#'
#' @param eem An EEM object.
#' @param pix_bins number of adjacent emission increments (pixels on a CCD fluoromter) to group together.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr group_split
#'
#' @export
#'
eem_bin_emission <- function(eem, pix_bins = 2){
  ## testvars
  # eem <- eem_test
  # pix_bins = 2
  # average = TRUE
  ## Catch and stop if the number of increments can't be evenly summed as pixels
  if(!length(eem$em)/pix_bins %% 2 != 0){
    stop(paste0("cannot bin emission increments as the number of increments (",length(eem$em),") is not evenly divisible by the number of pixel bins (",pix_bins,")."))
  }
  ## Pull out individual emission scans as a list, then bin them.
  eem_df_ls <- eem %>%
    as.data.frame(., gather = TRUE) %>%
    group_by(ex) %>%
    group_split(.) %>%
    lapply(., bin_emscan, pix_bins = pix_bins)
  eem_new <- eemdf_to_eem(eem_df_ls %>% rlist::list.rbind(),
                          file = eem$file,
                          sample = eem$sample,
                          location = eem$location,
                          gathered = TRUE)
  return(eem_new)
}
