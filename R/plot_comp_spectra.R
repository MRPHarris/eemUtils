#' Plot PARAFAC model Ex and Em mode (B and C) loadings
#'
#' @description Basic lineplot of Ex and Em spectra extracted from PARAFAC model components.
#'
#' @param pfmodel A PARAFAC model object.
#' @param eemlist An eemlist object.
#' @param component Integer. Which PARAFAC component to plot?
#' @param denorm TRUE/FALSE denormalise the loadings based upon eemlist (not parafac) fmax values?
#' @param labels TRUE/FALSE extract numeric values from sample names and use these as labels. values/10 by default to match a lexicographic scheme.
#' @param label_threshold Integer. Values lying above this number multiplied by the mean will be given labels via geom_text
#'
#' @export
#'
plot_comp_spectra <- function(pfmodel, comp){
  em_spec <- pfmodel$B %>%
    data.frame() %>%
    select(all_of(paste0('Comp.',comp))) %>%
    rownames_to_column('wavelength') %>%
    mutate(name = 'em') %>%
    rename_with(.cols = 2, ~"value") %>%
    mutate_at(vars(wavelength, value), as.numeric) %>%
    mutate(value = (value/max(value)))
  ex_spec <- pfmodel$C %>%
    data.frame() %>%
    select(all_of(paste0('Comp.',comp))) %>%
    rownames_to_column('wavelength') %>%
    mutate(name = 'ex') %>%
    rename_with(.cols = 2, ~"value") %>%
    mutate_at(vars(wavelength, value), as.numeric) %>%
    mutate(value = (value/max(value)))
  spec <- rbind(em_spec,ex_spec)
  p <- ggplot() +
    geom_line(data = spec, aes(x = wavelength, y = value, group = name)) +
    scale_x_continuous(expand = c(0,0)) +
    theme_cowplot(12)
  p
}
