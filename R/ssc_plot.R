
#' Plot two spectra and show congruence metric data
#'
#' @description plots two matrices (e.g. spectra) in a similar fashion to figure 2 from Wunsch et al., 2019.
#'
#' @param mat1 a matrix
#' @param mat2 a matrix
#' @param normalise TRUE/FALSE to normalise the matrices before comparison and subsequent plotting.
#' @param labels Optional character vectors to label the spectra with. Correspond to mat1, mat2.
#'
#' @export
#'
ssc_plot <- function(mat1, mat2, normalise = TRUE, labels = NULL){
  if(!is.null(labels)){
    mat_labels <- labels
  } else {
    mat_labels <- c(deparse(substitute(mat1)), deparse(substitute(mat2)))
  }
  if(!all(names(mat1) == names(mat2))){
    message("Warning: there is a misalignment in the matrix wavelengths.")
  }
  if(isTRUE(normalise)){
    mat1 <- mat1/max(mat1)
    mat2 <- mat2/max(mat2)
  }
  # Get SSC and associated vars
  ssc_fp <- ssc_more(mat1 = mat1, mat2 = mat2, tcc = FALSE)
  tcc_fp <- ssc(mat1 = mat1, mat2 = mat2, tcc = TRUE)
  # format data for plotting.
  SSC_table <- data.frame(matrix(cbind(mat1,mat2,rownames(mat1), rownames(mat2)), ncol = 4))
  colnames(SSC_table) <- c("mat1_value","mat2_value","mat1_wavelength", "mat2_wavelength")
  SSC_table <- SSC_table %>%
    pivot_longer(everything(), names_sep = "_", names_to = c("mat",".value")) %>%
    mutate_at(vars(value,wavelength), as.numeric)
  # colours
  clrs_fill <- c("#6b6b6b","#3e80a8")
  clrs_line <- c("#383838","#234961")
  # plotting
  plt <- ggplot() +
    geom_line(data = SSC_table, aes(x = wavelength, y = value, group = mat, colour = mat, linetype = mat), size = 0.8) +
    theme_cowplot(12) +
    geom_area(data = SSC_table, aes(x = wavelength, y = value, fill = mat, group = mat), alpha = 0.5, position = "identity") +
    scale_linetype_manual(values = c(1, 2), labels = mat_labels) +
    scale_color_manual(values = clrs_line, labels = mat_labels) +
    scale_fill_manual(values = clrs_fill, labels = mat_labels) +
    scale_y_continuous(expand = c(0,0), limits = c(0,(max(SSC_table$value)+(max(SSC_table$value)/10)))) +
    scale_x_continuous(expand = c(0,0)) +
    labs(fill = "Spectra", color = "Spectra", linetype = "Spectra",
         x = "Wavelength") +
    annotate(geom = "text", x = (max(SSC_table$wavelength)-(max(SSC_table$wavelength)/10)), y = (max(SSC_table$value)-(max(SSC_table$value)/10)),
             label = paste0("TCC = ",    format(round(tcc_fp,3), nsmall = 3), "\n",
                            "SSC = ",    format(round(ssc_fp[1,],3), nsmall = 3), "\n",
                            "\u03B1 = ", format(round(ssc_fp[2,],3), nsmall = 3),"\n",
                            "\u03B2 = ", format(round(ssc_fp[3,],3), nsmall = 3)),
             hjust = 1)
  if(isTRUE(normalise)){
    plt <- plt +
      labs(y = "Normalised value")
  }
  plt
}
