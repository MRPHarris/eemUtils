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
