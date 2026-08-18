#' Multiply the loadings values of PARAFAC components by normalisation factors
#'
#' @description A PARAFAC model fed with normalised eems outputs by default loadings that are less useful for direct
#'       sample-to-sample comparison. This can be remedied by multiplying the loadings by each eem fmax value (i.e. the
#'       normalisation factors used to normalise the eems originally).
#'
#' @param pfmodel a PARAFAC model object, typically an individual ouptut from staRdom::eem_parafaC()
#' @param eemlist a list of eem objects compliant with the staRdom/eemR framework
#' @param type short or long format. Short by default. Long format data works better for grouping in ggplot
#'
#' @export
#'
extrpf_loadings_denorm <- function(pfmodel, eemlist, type = "short"){
  maxvals <- eemlist_fmax_values(eemlist)
  fm <- extrpf_loadings(pfmodel)
  loadings_frame <- fm[,2:ncol(fm)] %>% as.data.frame()
  newframe <- apply(loadings_frame, 2, function(col) {
    col * maxvals
  }) %>% data.frame()
  newframe <- newframe %>% `colnames<-`(c(paste0("Comp.",
                                                 seq(1, ncol(newframe), 1)))) %>% `rownames<-`(unlist(lapply(eemlist,
                                                                                                             "[[", "sample"))) %>% rownames_to_column(var = "sample")
  if (type == "short") {
    newframe
  }
  else if (type == "long") {
    newframe <- newframe %>% pivot_longer(cols = c(2:ncol(newframe)))
  }
  else {
    stop("Unknown 'type'. Please input either 'short' or 'long'")
  }
}
