#' Plot spectra from an EEM at a supplied Ex/Em position
#'
#' @description A simple ggplot2 lineplot producer, using slice_eem to extract Ex and Em
#'      spectra at the given coordinates.
#'
#' @param eem a single EEM object, compliant with the staRdom/eemR framework
#' @param ex numeric excitation wavelength value
#' @param em numeric emission wavelength value
#'
#' @export
#'
plot_eem_spectra <- function(eem, ex, em){
  # var checks
  if(class(eem) != 'eem'){
    stop("Please pass the function an object of class 'eem'")
  }
  spectra <- slice_eem(eem = eem, ex = ex, em = em)
  p <- ggplot() +
    geom_line(data = spectra, aes(x = wavelength, y = intensity, group = name, linetype = name)) +
    # geom_area(data = spectra, aes(x = wavelength, y = intensity, fill = name), alpha = 0.5) +
    labs(x = c("Wavelength (nm)"),
         y = "Value") +
    scale_linetype_manual(values = c('solid','dashed')) +
    theme_cowplot(12) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme(legend.position = 'none')
  p
}
