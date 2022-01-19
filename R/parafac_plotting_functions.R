# Functions for the plotting of PARAFAC-related datasets.

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

#' Plot PARAFAC component modeled EEM in 3D, with a reversed excitation axis.
#'
#' @description A tweaked staRdom::eempf_comps3D() function that accounts for
#'      the excitation axis being incorrectly plotted in the original. Seems to
#'      be a common feature for Aqualog .DAT or .csv originating data. Not quite
#'      sure where the flipped axis enters the equation!
#'
#' @param pfmodel A PARAFAC model object.
#' @param which NULL or numeric for which component from the PARAFAC model to plot.
#'
#' @export
#'
eempf_comps3D_revex <- function (pfmodel, which = NULL){
  data <- pfmodel %>% eempf_comp_mat()
  z <- lapply(data, function(mat) {
    mat %>% data.frame() %>% spread(em, value) %>% remove_rownames() %>%
      column_to_rownames("ex") %>% as.matrix()
  })
  ex <- lapply(data, function(mat) {
    mat$ex %>% rev() %>% unique() %>% as.numeric() %>% as.vector()
  })
  em <- lapply(data, function(mat) {
    mat$em %>% unique() %>% as.numeric() %>% as.vector()
  })
  scene <- list(xaxis = list(title = "em"), yaxis = list(title = "ex"))
  lapply(1:length(ex), function(comp) {
    if (is.null(which) | comp %in% which) {
      plotly::plot_ly(x = em[[comp]], y = ex[[comp]], z = z[[comp]],
                      colors = rainbow(12)[9:1]) %>% plotly::layout(scene = scene) %>%
        plotly::add_surface()
    }
  })
}

#' Plot PARAFAC model loadings versus per-sample Fmax at the component peak position.
#'
#' @description A simple scatter plot of PARAFAC component loadings versus Fmax data. Optional labelling,
#'      derived through extraction of numeric values in the sample names.
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
fmax_corrplot <- function(pfmodel, eemlist, component = 1, denorm = TRUE, labels = FALSE, label_threshold = 2){
  if(isTRUE(denorm)){
    loading <- extrpf_loadings_denorm(pfmodel, eemlist) %>%
      dplyr::select(paste0("Comp.",component),sample) %>%
      'colnames<-'(c(paste0("loading"),"sample"))
  } else {
    loading <- extrpf_loadings(pfmodel) %>%
      dplyr::select(paste0("Comp.",component),sample) %>%
      'colnames<-'(c(paste0("loading"),"sample"))
  }
  fmax <- extrpf_fmax(pfmodel, eemlist) %>%
    dplyr::select(paste0("Comp.",component),sample) %>%
    'colnames<-'(c(paste0("fmax"),"sample"))
  # Bind, pivot
  frame <- merge(loading,fmax)
  if(isTRUE(labels)){
    frame$index <- as.numeric(regmatches(frame$sample,gregexpr("[[:digit:]]+", frame$sample)))
  }
  # lm equation
  df <- frame[,2:3]
  colnames(df) <- c("x","y")
  my.formula <- y ~ x
  # ggplot
  plt <- ggplot(data = frame, aes(x = loading, y = fmax)) +
    geom_smooth(method = 'lm') +
    stat_poly_eq(data = df, formula = my.formula,
                 aes(x = x, y = y, label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                 parse = TRUE) +
    geom_point(shape = 21, color = "black", fill = "white", size = 1.5) +
    scale_x_continuous(expand = c(0,0), limits = c(0-max(frame$loading/50), ceiling(max(frame$loading)), breaks = seq(0,ceiling(max(frame$loading)),1))) +
    scale_y_continuous(expand = c(0,0), limits = c(0-max(frame$fmax/50),ceiling(max(frame$fmax))), breaks = seq(0,ceiling(max(frame$fmax)),1)) +
    theme_cowplot(12)
  if(isTRUE(labels)){
    # label samples only greater than certain value
    if(label_threshold == 0){
      plt <- plt +
        geom_text(data = frame, aes(x = loading, y = fmax, label = ifelse(loading > 0 | fmax > 0 ,as.character(index/10),'')),
                  size = 2, hjust = -0.3, vjust = -0.3)
    } else {
      plt <- plt +
        geom_text(data = frame, aes(x = loading, y = fmax, label = ifelse(loading > (label_threshold*mean(loading)) | fmax > (label_threshold*mean(fmax)) ,as.character(index/10),'')),
                  size = 3, hjust = -0.3, vjust = -0.3)
    }
  }
  plt
}


