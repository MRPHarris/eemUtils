#' Plot an EEM object in 3D using plotly.
#'
#' @description A conversion of staRdom::eempf_comps3D(), applied to EEM objects.
#'
#' @param eem An eem object, compliant with the eemR/staRdom packages.
#' @param origin How was this eem imported? One of 'csv' (generic) or 'DAT' (Aqualog import)
#'
#' @export
#'

plot_eem_3D <- function (eem, origin = "csv"){
  if(!(origin == "DAT" || origin == "csv")){
    warning("please specify an origin for the eem: either 'csv' or 'DAT'")
  }
  if(class(eem) == "eem"){
    eem2df_nested <- function (eem, gather = TRUE){
      eem_df <- eem[["x"]] %>% data.frame()
      colnames(eem_df) <- eem[["ex"]]
      rownames(eem_df) <- eem[["em"]]
      if (gather == TRUE){
        eem_df <- eem_df %>% tibble::rownames_to_column("em") %>%
          gather(ex, value, -em)}
      eem_df
    }
    data <- eem %>% eem2df_nested()
  } else(warning("please use an object of class 'eem'"))
  z <- data %>% data.frame() %>% spread(em, value) %>% remove_rownames() %>% column_to_rownames("ex") %>% as.matrix()
  ex <- data$ex %>% unique() %>% as.numeric() %>% as.vector()
  if(origin == "csv"){
    ex <- rev(ex)
  }
  em <- data$em %>% unique() %>% as.numeric() %>% as.vector()
  scene <- list(
    xaxis = list(
      title = "em"),
    yaxis = list(
      title = "ex",
      autorange = "reversed"))
  plotly::plot_ly(x = em, y = ex, z = z, colors = rainbow (12)[9:1]) %>%
    plotly::layout(scene = scene) %>%
    plotly::add_surface()
}
