#' Extract the modelled EEMs from a given PARAFAC model as an eemlist.
#'
#' @description Pulls out the modelled EEMs from a PARAFAC model, in the standard eemlist format,
#'      comprised of eem objects compliant with the staRdom/EEM/eemR framework.
#'
#' @param pfmodel a PARAFAC model object containing one or more components
#'
#' @export
#'
extrpf_eems <- function(pfmodel){
  # Code based upon staRdom::eempft_comp_mat()
  eemlist <- lapply(seq(1:ncol(pfmodel$A)), function(comp){
    m <- matrix(pfmodel$B[, comp]) %*% t(matrix(pfmodel$C[,comp])) %>% data.frame()
    colnames(m) <- pfmodel$C %>% rownames()
    rownames(m) <- pfmodel$B %>% rownames()
    m_eem <- eemUtils::eemdf_to_eem(m, file = NULL, sample = colnames(pfmodel$B)[comp], location = NULL, gathered = FALSE)
    m_eem
  })
  names(eemlist) <- colnames(pfmodel$B)
  eemlist <- eemlist %>% 'class<-'(c('eemlist'))
  eemlist
}
