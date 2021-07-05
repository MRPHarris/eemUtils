# These functions relate to the plotting of excitation-emission matrices.

#' Visually check individual EEMs within an eemlist.
#'
#' @description Display and/or export a single EEM, for inspection/checking purposes.
#'
#' @param eemlist A list of EEMs, compliant with the eemR/staRdom framework.
#' @param contour TRUE/FALSE to show contours on plot.
#' @param eem_number The number of the eem to plot.
#' @param output_dir Directory to send a .png of the plot to. If no export is desired, leave as NULL.
#'
#' @export
#'
check_eem <- function(eemlist, contour = TRUE, eem_number = 1, output_dir = NULL){
  EEM_Name = eemlist[[eem_number]][["sample"]]
  print(EEM_Name)
  message(paste0("EEM ",eem_number,"/",length(eemlist)))
  Target_eemlist <- vector("list", 1)
  Target_eemlist[[1]] <- eemlist[[eem_number]]
  output_dir = output_dir
  if(!(is.null(output_dir))){
    png(paste0(output_dir,EEM_Name,".png"), units="cm", width=18, height=14, res=300)
    print(eem_overview_plot(Target_eemlist, contour = contour)) # change here to output contours in .png EEMs
    dev.off()
  } else(
    message("no output directory defined; png not exported")
  )
  print(eem_overview_plot(Target_eemlist, contour = contour))
  Target_eemlist <<- Target_eemlist
}

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
