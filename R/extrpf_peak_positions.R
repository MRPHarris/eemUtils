#' Return the ex/em position of the peak maxima for a PARAFAC modeled component.
#'
#' @description Uses extrpf_peak_spectra to obtain the single ex/em coordinate pair for the maxima of the
#'     specified PARAFAC component.
#'
#' @param pfmodel A PARAFAC model object.
#' @param component NULL for all components, or the numeric identifier of a single component.
#'
#' @export
#'
extrpf_peak_positions <- function(pfmodel, component = NULL){
  if(!is.null(component) || is.numeric(component)){
    stop("'component' must be NULL or numeric")
  }
  peak_spectra = extrpf_peak_spectra(pfmodel = pfmodel,
                                     component = component)
  peakpositions = peak_spectra %>%
    dplyr::select(max_ex,max_em) %>%
    distinct()
  peakpositions
}
