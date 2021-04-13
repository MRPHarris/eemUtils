#' Plot PARAFAC component modeled EEMs.
#'
#' @description Plot the modeled EEMs from one or more PARAFAC models. Takes an
#'      output from extrpf_spectra_or_eems(type = 1).
#'
#' @param modeled_EEMs An output from extrpf_spectra_or_eems(type = 1).
#' @param contour TRUE/FALSE to display contours on the EEM plot/s.
#'
#' @export
#'

# plot modeled eems. Takes an output from extrpf_spectra_or_eems(type = 1).
plot_extrpf_eems <- function(modeled_EEMs, contour = FALSE){
  tab <- modeled_EEMs
  colpal <- rainbow(75)[53:1]
  fill_max <- tab$value %>% max(na.rm = TRUE)
  vals <- seq(from = 0, to = fill_max, length.out = length(colpal))
  vals <- (vals - min(vals))/diff(range(vals))
  colpal <- rainbow(75)[53:1]
  # Diffs?
  diffs <- tab %>%
    select(-value, -comps) %>%
    gather("spec", "wl", -comp, -modname) %>%
    group_by(comp, modname, spec) %>%
    unique() %>%
    summarise(slits = diff(wl) %>% n_distinct()) %>%
    .$slits != 1
  # Plotting
  plot <- tab %>%
    ggplot(aes(x = ex, y = em, z = value))
  # Standard tiles: diffs present
  if (any(diffs)) {
    plot <- plot +
      layer(mapping = aes(colour = value, fill = value),
            geom = "tile",
            stat = "identity",
            position = "identity")
  }
  # Rasters: diffs not present
  else {
    plot <- plot +
      layer(mapping = aes(fill = value),
            geom = "raster",
            stat = "identity",
            position = "identity")
  }
  plot <- plot +
    scale_fill_gradientn(colours = colpal,
                         values = vals,
                         limits = c(tab$value %>% min(na.rm = TRUE), fill_max),
                         aesthetics = c("fill", "colour")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = "Excitation (nm)", y = "Emission (nm)") +
    facet_grid(comp ~ modname)
  if (isTRUE(contour)) {
    plot <- plot +
      geom_contour(colour = "black",size = 0.3)
  }
  plot
}
