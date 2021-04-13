#' Extract fluorescence intensity loadings from a list of PARAFAC models.
#'
#' @description Extract the per-sample modeled intensity loadings for each component
#'       from a list of PARAFAC models. Optional numeric indexing of sample names
#'       using regmatching, assuming there are numerals to extract.
#'
#' @param pfmodel_list A list of PARAFAC models, typically an output from stardom::eem_parafac()
#' @param by_index TRUE/FALSE index the loadings by the sample names via registry matching. Requires numerals in the sample names.
#'
#' @export
#'

extrpf_loadings <- function(pfmodel_list, by_index = FALSE){
  loadings_list <- vector(mode = "list", length = length(pfmodel_list))
  pfmodel_name <- deparse(substitute(pfmodel_list))
  for (i in seq_along(pfmodel_list)){
    model <- pfmodel_list[[i]][["A"]]
    LoadingsA <- as.data.frame(model)
    Loadings <- as.data.frame(model)
    comp_num <- vector(mode = "list", length = ncol(Loadings))
    for (i in seq_along(comp_num)){
      cstat <- data.frame(matrix((paste0("C",i)),nrow(Loadings),1))
      Loadings <- cbind.data.frame(Loadings,cstat)
      col_index = i+(ncol(LoadingsA))
      colnames(Loadings)[col_index] <- paste0("C",i,"_Status")
    }
    data.table::setDT(Loadings, keep.rownames = TRUE)[]
    colnames(Loadings)[1] <- c("Sample_Name")
    if(isTRUE(by_index)){
      index <- as.numeric(regmatches(Loadings$Sample_Name, gregexpr("[[:digit:]]+", Loadings$Sample_Name)))
      Loadings <- cbind(Loadings, index)
    }
    loadings_list[[i]] <- Loadings
  }
  loadings_list <<- loadings_list
  assign(paste0("Loadings_",pfmodel_name), loadings_list, envir = parent.frame())
  rm(loadings_list, envir = parent.frame())
  message("Success! Loadings data extracted as Loadings_",pfmodel_name,".")
}
