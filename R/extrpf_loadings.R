#' Extract fluorescence intensity loadings from a PARAFAC model.
#'
#' @description Extract the per-sample modeled intensity loadings for each component
#'       from a PARAFAC models. Optional numeric indexing of sample names
#'       using regmatching, assuming there are numerals to extract. Use lapply for
#'       multiple models.
#'
#' @param pfmodel A PARAFAC model object containing any number of components
#' @param by_index TRUE/FALSE index the loadings by the sample names via registry matching. Requires numerals in the sample names
#'
#' @export
#'
extrpf_loadings <- function(pfmodel, by_index = FALSE, fill_stat = FALSE){
  model <- pfmodel[["A"]]
  LoadingsA <- as.data.frame(model)
  Loadings <- as.data.frame(model)
  comp_num <- vector(mode = "list", length = ncol(Loadings))
  if(isTRUE(fill_stat)){
    for (i in seq_along(comp_num)) {
      cstat <- data.frame(matrix((paste0("C", i)),nrow(Loadings), 1))
      Loadings <- cbind.data.frame(Loadings, cstat)
      col_index = i + (ncol(LoadingsA))
      colnames(Loadings)[col_index] <- paste0("C",i, "_Status")
    }
  }
  # data.table::setDT(Loadings, keep.rownames = TRUE)[]
  Loadings <- Loadings %>%
    rownames_to_column("sample") %>%
    dplyr::select(sample, everything())
  colnames(Loadings)[1] <- c("sample")
  if (isTRUE(by_index)) {
    index <- as.numeric(regmatches(Loadings$sample,gregexpr("[[:digit:]]+", Loadings$sample)))
    Loadings <- cbind(Loadings, index)
  }
  Loadings
}
