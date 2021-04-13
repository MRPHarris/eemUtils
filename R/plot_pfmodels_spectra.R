#' Plot the component spectra intersecting at component maxima from a list of PARAFAC models.
#'
#' @description This is an unpacked staRdom::eempf_plot_comps(type = 2).
#'
#' @param pfmodel_list PARAFAC model list.
#' @param eemlist list of EEMs used to produce the pfmodels.
#' @param title A title to add to the plot.
#'
#' @export
#'

plot_pfmodels_spectra <- function (pfmodel_list, eemlist, title){
  pfres <- pfmodel_list
  names = TRUE
  c <- pfres %>% lapply(eempf_comp_mat)
  colpal <- rainbow(75)[53:1]
  if (is.null(names(c))){
    names(c) <- paste0("model", seq(1:length(c)))
  }
  # Main results table, pre main rearrange
  tab <- lapply(1:length(c), function(n) {
    c1 <- c[[n]]
    mod_name <- names(c)[n]
    nc1 <- length(c1)
    nc2 <- 0
    lapply(c1, function(c2) {
      nc2 <<- nc2 + 1
      c2 <- c2 %>%
        mutate(comps = nc1,
               comp = paste0("Comp.", nc2),
               modname = mod_name)
    }) %>%
      bind_rows()
  }) %>%
    bind_rows() %>%
    mutate(modname = factor(modname, levels = names(c))) %>%
    mutate_at(vars(ex, em, value), as.numeric)
  fill_max <- tab$value %>% max(na.rm = TRUE)
  vals <- seq(from = 0, to = fill_max, length.out = length(colpal))
  vals <- (vals - min(vals))/diff(range(vals))
  # Creating object for spectra plot
  max_nm <- as.numeric(max(eemlist[[1]]$em, na.rm = TRUE))
  plot_data <- tab %>%
    group_by(modname, comp) %>%
    mutate(max_pos = which.max(value), max_em = em[max_pos], max_ex = ex[max_pos]) %>%
    mutate(exn = ifelse(em == max_em, ex, NA), emn = ifelse(ex == max_ex, em, NA)) %>%
    filter(!is.na(emn) | !is.na(exn)) %>%
    ungroup()
  ggplot() +
    geom_line(data = plot_data, aes(x = exn, y = value), colour = "#000000", group = "excitation", na.rm = TRUE) + # line for first component
    geom_area(data = plot_data, aes(x = exn, y = value, fill = "excitation"), alpha = 0.6, group = "excitation") + # area fill for first component
    geom_line(data = plot_data, aes(x = emn, y = value), colour = "#000000", group = "emission", na.rm = TRUE) + # line for second component
    geom_area(data = plot_data, aes(x = emn, y = value, fill = "emission"), alpha = 0.6, group = "emission") +
    labs(
      title = title,
      x = "Wavelength (nm)",
      y = "Loading") +
    scale_fill_manual(
      values = c("darkblue","steelblue"), # CHANGE FILL COLOURS HERE
      #values = c("#033d04","darkgreen"),
      #values = ExtendedHotCold(num_components_val),
      name = "Component Spectra") +
    scale_x_continuous(
      expand = c(0,0),
      limits = c(235,max_nm),
      breaks = seq(250,max_nm,50),
      position = "bottom") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.line = element_line(),
      plot.title = element_text(hjust = 0.5),
      legend.title = element_text(size=10, face = "bold"),
      legend.position = c(0.1,0.15),
      legend.background = element_rect(fill = "white", size = 0.5, linetype = "solid", colour = "white"),
      panel.background = element_rect(fill = NA, color = "#d4d4d6"),
      panel.grid.major.x = element_line(size = 0.3, linetype = 'dashed', colour = '#d4d4d6')) +
    facet_grid(
      comp ~ modname)
}
