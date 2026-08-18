#' Extract the per-EEM percentage contributions of each component
#'
#' @description Return the percent contribution of each component to the modelled fluorescence response of each sample, using
#'      the modelled loadings
#'
#' @param pfmodel a PARAFAC model object
#' @param eemlist supply the eemlist used for denormalisation. Only necessary if denormalise is set to TRUE.
#' @param denormalise TRUE/FALSE to denormalise the loadings prior to the percent calculation
#'
#' @export
#'
extrpf_loadings_percent <- function(pfmodel, eemlist, denormalise = FALSE){
  if(isTRUE(denormalise)){
    loadings <- extrpf_loadings_denorm(pfmodel, eemlist)
  } else {
    loadings <- extrpf_loadings(pfmodel)
  }
  FI_totals <- rowSums(loadings[,2:ncol(loadings)])
  pct_contrib <- loadings[,2:ncol(loadings)] %>%
    mutate_all(.,function(col){(col/FI_totals)*100})
  pct_contrib$sample <- loadings$sample
  pct_contrib <- pct_contrib %>%
    dplyr::select('sample',everything())
  pct_contrib
}
