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

eempf_comps3D_revex <- function (pfmodel, which = NULL)
{
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
