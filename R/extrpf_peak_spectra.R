#' Extract spectra at the peak wavelength position.
#'
#' @description Extract the spectra at the peak for a given component.
#'
#' @param pfmodel A PARAFAC model object
#' @param component NULL or numeric for specific components.
#'
#' @export
#'
extrpf_peak_spectra <- function(pfmodel, component = 1){
  pfres <- pfmodel # Code adapted from staRdom package; naming will reflect this.
  names = TRUE
  c <- eempf_comp_mat(pfres) # get matrices of all the comps.
  if (is.null(names(c))) {
    names(c) <- paste0("model", seq(1:length(c)))
  }
  # Performing for only one model.
  c1 <- c # holdover code; lapply wrap removed here.
  nc1 <- length(c1)
  nc2 <- 0
  # for each component, pull out the data
  tab <- lapply(c1, function(c2) {
    nc2 <<- nc2 + 1
    c2 <- c2 %>% mutate(comps = nc1, comp = paste0(nc2))
  }) %>% bind_rows() %>%
    mutate_at(vars(ex, em, value, comp), as.numeric) # numeric vars
  # collate, organise, extract max spectra.
  comp_spectra <- tab %>%
    group_by(comp) %>%
    mutate(max_pos = which.max(value), max_em = em[max_pos], max_ex = ex[max_pos]) %>%
    mutate(exn = ifelse(em == max_em, ex, NA), emn = ifelse(ex == max_ex, em, NA)) %>%
    dplyr::filter(!is.na(emn) | !is.na(exn)) %>%
    ungroup() %>%
    mutate_at(vars(ex, em, value, comp), as.numeric)# %>%
  # dplyr::filter(comp == component)
  if(!is.null(component)){
    comp_spectra <- comp_spectra %>%
      dplyr::filter(comp == component)
  }
  comp_spectra
}
