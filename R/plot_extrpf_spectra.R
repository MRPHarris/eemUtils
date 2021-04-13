#' Plot PARAFAC component peak maxima spectra.
#'
#' @description Plot the modeled spectra from one or more PARAFAC models. Takes
#'      an output from extrpf_spectra_or_eems(type = 2).
#'
#' @param modeled_spectra An output from extrpf_spectra_or_eems(type = 2).
#' @param eemlist The list of EEMs used to produce the PARAFAC.
#' @param title A character vector, to be displayed as a title on the plot.
#'
#' @export
#'

plot_extrpf_spectra <- function(modeled_spectra, eemlist, title){
  spectra <- modeled_spectra #shorten input
  Component_names <- list(
    'Comp.1'="Component 1",
    'Comp.2'="Component 2",
    'Comp.3'="Component 3",
    'Comp.4'="Component 4",
    'Comp.5'="Component 5",
    'Comp.6'="Component 6",
    'Comp.7'="Component 7",
    'Comp.8'="Component 8",
    'Comp.9'="Component 9")
  Mod_names <- list(
    'model1'="Model 1",
    'model2'="Model 2",
    'model3'="Model 3",
    'model4'="Model 4",
    'model5'="Model 5",
    'model6'="Model 6",
    'model7'="Model 7",
    'model8'="Model 8",
    'model9'="Model 9")
  compmod_labeller <- function(variable,value){
    if (variable=='comp') {
      return(Component_names[value])
    } else {
      return(Mod_names[value])
    }
  }
  #num_components_val = 2
  #ExtendedHotCold = colorRampPalette(brewer.pal(11, "RdYlBu"))  # Optional adaptive ramped palette
  max_nm <- as.numeric(max(eemlist[[1]]$em, na.rm = TRUE))
  plot <- ggplot() # init ggplot object
  plot_exp <- plot +
    geom_line(data = spectra, aes(x = exn, y = value), colour = "#000000", group = "excitation", na.rm = TRUE) + # line for first component
    geom_area(data = spectra ,aes(x = exn, y = value, fill = "excitation"), alpha = 0.5, group = "excitation") + # area fill for first component
    geom_line(data = spectra, aes(x = emn, y = value), colour = "#000000", group = "emission", na.rm = TRUE) + # line for second component
    geom_area(data = spectra, aes(x = emn, y = value, fill = "emission"), alpha = 0.5, group = "emission") +
    labs(
      title = title,
      x = "Wavelength (nm)",
      y = "Loading") +
    scale_fill_manual(
      values = c("darkblue","steelblue"), # CHANGE FILL COLOURS HERE
      #values = ExtendedHotCold(num_components_val),
      name = "Component Spectra") +
    scale_x_continuous(
      expand = c(0,0),
      limits = c(235,max_nm),
      breaks = seq(250,max_nm,50),
      position = "bottom") +
    scale_y_continuous(
      expand = c(0,0)) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.line = element_line(),
      plot.title = element_text(hjust = 0.5),
      legend.title = element_text(size=10, face = "bold"),
      legend.position = c(0.8,0.4),
      legend.background = element_rect(fill = "white", size = 0.5, linetype = "solid", colour = "white"),
      panel.background = element_rect(fill = NA, color = "#d4d4d6"),
      panel.grid.major.x = element_line(size = 0.3, linetype = 'dashed', colour = '#d4d4d6')) +
    facet_grid(
      comp ~ modname,
      labeller = compmod_labeller)
  plot_exp <<- plot_exp
  print(plot_exp)
}
