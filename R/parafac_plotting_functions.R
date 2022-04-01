# Functions for the plotting of PARAFAC-related datasets.

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

#' Plot PARAFAC model loadings versus per-sample Fmax at the component peak position.
#'
#' @description A simple scatter plot of PARAFAC component loadings versus Fmax data. Optional labelling,
#'      derived through extraction of numeric values in the sample names. Given fmax is just rescaled loadings data,
#'      this will typically be perfectly correlated. So, the recommended setting for 'type' is 'peakpick', in which case
#'      extrpf_fmax() will fetch the intensity at the component maxima in each EEM. The resulting data will be a measure of
#'      sample-to-model homogeneity, reflecting intra-component variance, underfitting, or spectral overlap at the peak coordinates.
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
fmax_peakpick_corrplot <- function(pfmodel, eemlist, component = 1, denorm = TRUE, labels = FALSE, label_threshold = 2){
  fmax <- extrpf_fmax(pfmodel, eemlist, component = component, type = 'fmax', denormalise = denorm) %>%
    dplyr::select(paste0("Comp.",component),sample) %>%
    'colnames<-'(c(paste0("fmax"),"sample"))
  peakint <- extrpf_fmax(pfmodel, eemlist, component = component, type = 'peakpick') %>%
    dplyr::select(paste0("Comp.",component),sample) %>%
    'colnames<-'(c(paste0("peakint"),"sample"))
  # Bind, pivot
  frame <- merge(fmax,peakint)
  if(isTRUE(labels)){
    frame$index <- as.numeric(regmatches(frame$sample,gregexpr("[[:digit:]]+", frame$sample)))
  }
  # lm equation
  df <- frame[,2:3]
  colnames(df) <- c("x","y")
  my.formula <- y ~ x
  # ggplot
  plt <- ggplot(data = frame, aes(x = fmax, y = peakint)) +
    geom_smooth(method = 'lm') +
    stat_poly_eq(data = df, formula = my.formula,
                 aes(x = x, y = y, label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                 parse = TRUE) +
    geom_point(shape = 21, color = "black", fill = "white", size = 1.5) +
    scale_x_continuous(expand = c(0,0), limits = c(0-max(frame$fmax/50), ceiling(max(frame$fmax)), breaks = seq(0,ceiling(max(frame$fmax)),1))) +
    scale_y_continuous(expand = c(0,0), limits = c(0-max(frame$peakint/50),ceiling(max(frame$peakint))), breaks = seq(0,ceiling(max(frame$peakint)),1)) +
    theme_cowplot(12) +
    labs(x = "Fmax", y = "Component peak intensity")
  if(isTRUE(labels)){
    # label samples only greater than certain value
    if(label_threshold == 0){
      plt <- plt +
        geom_text(data = frame, aes(x = fmax, y = peakint, label = ifelse(loading > 0 | fmax > 0 ,as.character(index/10),'')),
                  size = 2, hjust = -0.3, vjust = -0.3)
    } else {
      plt <- plt +
        geom_text(data = frame, aes(x = fmax, y = peakint, label = ifelse(fmax > (label_threshold*mean(fmax)) | peakint > (label_threshold*mean(peakint)) ,as.character(index/10),'')),
                  size = 3, hjust = -0.3, vjust = -0.3)
    }
  }
  plt
}

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

