#' Extract per-sample residual EEM objects
#'
#' Extract the 'leftover' fluorescence (i.e. residuals) in each EEM passed to a given PARAFAC model. Areas of
#'     high (negative) fluorescence indicate regions where the PARAFAC model is under (over) fitted.
#'
#' @param pfmodel a single PARAFAC model object returned from staRdom::eem_parafac containing any number of components.
#' @param eemlist The eemlist object passed to the above PARAFAC model.
#'
#' @importFrom dplyr group_by
#' @importFrom dplyr filter
#' @importFrom dplyr select
#'
#' @export
#'
extrpf_residual_eems <- function(pfmodel,
                                 eemlist){
  resids <- extrpf_residuals(pfmodel = pfmodel,
                             eem_list = eemlist)
  res_frame <- resids %>%
    dplyr::filter(type == 'residual') %>%
    dplyr::select(-type) %>%
    group_by(sample)
  res_samp_list <- split(x = res_frame, f = res_frame$sample)
  # Loop because lapply is finnicky with indices
  res_samp_eems <- vector('list',length = length(eemlist)) %>% 'class<-'(c('eemlist'))
  for(e in seq_along(res_samp_eems)){
    eem_it <- res_samp_list[[e]]
    eemdf_it <- eemdf_to_eem(eemdf = eem_it,
                             file = NULL,
                             sample = paste0(unlist(lapply(eemlist,"[[",'sample'))[e],"_residuals"),
                             location = NULL,
                             gathered = TRUE)
    res_samp_eems[[e]] <- eemdf_it
  }
  return(res_samp_eems)
}
