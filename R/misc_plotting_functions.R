# Miscellaneous functions related to plotting that didn't match the other function compilation files

#' A simple textsize scaler for ggplot objects.
#'
#' @description Takes a set of text size alteration numbers and modifies a ggplot. Useful
#'      when attempting to get a ggplot object's text scaling right.
#'
#' @param plot a ggplot2 object.
#' @param multiplier numeric; a direct multiplier for specified text sizes.
#' @param default_title_size numeric; default text size for plot title.
#' @param default_axis_title_size numeric; default text size for axis titles.
#' @param default_axis_text_size numeric; default text size for axis text.
#' @param default_legend_title_size numeric; default text size for legend titles.
#' @param default_legend_text_size numeric; default text size for legend text
#'
#' @export
#'
ggplot_textsize_scale <- function(plot,
                                  multiplier = 0.8,
                                  default_title_size = 11,
                                  default_axis_title_size = 10,
                                  default_axis_text_size = 9,
                                  default_legend_title_size = 10,
                                  default_legend_text_size = 9){
  plot <- plot +
    theme(plot.title = element_text(size = default_title_size * multiplier),
          axis.text = element_text(size = default_axis_text_size * multiplier),
          axis.title = element_text(size = default_axis_title_size * multiplier),
          legend.text = element_text(size = default_legend_text_size * multiplier),
          legend.title = element_text(size = default_legend_title_size * multiplier))
  plot
}
