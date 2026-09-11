#' Bin adjacent emission increments/pixels.
#'
#' @description Average adjacent emission pixels to reduce noise in an emission spectra.
#'
#' @param em_spectra a long-form (gathered) emission spectra.
#' @param pix_bins number of adjacent emission increments (pixels on a CCD fluoromter) to group together.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr ungroup
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#'
#' @export
#'
bin_emscan <- function(em_spectra, pix_bins = 2, average = TRUE){
  # Bin spectra
  spectra_binned <- em_spectra %>%
    mutate(em = as.numeric(em)) %>%
    mutate(grp = ceiling(row_number() / pix_bins)) %>%
    group_by(grp) %>%
    summarise(value = sum(value)/pix_bins,
              ex = first(ex),
              em = mean(em),
              sample = first(sample)) %>%
    ungroup() %>% select(-grp)
  # Return
  return(spectra_binned)
}
