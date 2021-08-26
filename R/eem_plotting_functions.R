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
plot_eem_3D <- function (eem, origin = "DAT"){
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

#' Save all the EEMs within an eemlist object as .png images
#'
#' @description loop through an eemlist, saving EEM data as individual .png images.
#'
#' @param eemlist the target eemlist object compliant with the staRdom/eemR framework.
#' @param contour TRUE/FALSE to add contour to EEM.
#' @param output_dir full path to the destination folder. If NULL, images will be sent to the working directory.
#'
#' @export
#'
save_eemlist_pngs <- function(eemlist, contour = FALSE, output_dir = NULL){
  cyn = contour
  if(is.null(output_dir)){
    output_dir <- paste0(getwd(),"/")
  }
  for (i in seq_along(eemlist)){
    EEM_Name = eemlist[[i]][["sample"]]
    print(EEM_Name)
    message(paste0("EEM ",i,"/",length(eemlist)))
    Target_eemlist <- vector("list", 1)
    Target_eemlist[[1]] <- eemlist[[i]]
    png(paste0(output_dir,EEM_Name,".png"), units="cm", width=18, height=14, res=300)
    print(eem_overview_plot(Target_eemlist, contour = cyn)) # change here to output contours in .png EEMs
    dev.off()
  }
}

#' Extract a specific EEM from a set of eemlists
#'
#' @description when given a list of eemlists representing different stages in a processing
#'    framework, extract a single EEM from each eemlist in order to show the processing steps.
#'
#' @param list_of_eemlist a list vector, with each element comprising an eemlist compliant with the staRdom/eemR framework
#' @param which_eem integer - numerically index which eem within the eemlists?
#' @param output_dir NULL or full path to destination directory in which to save a .png image of the EEMs together.
#' @param contour TRUE/FALSE include contours in exported plots, assuming output_dir != NULL
#'
#' @export
#'
extract_procstep_eems <- function(list_of_eemlists, which_eem = 1, output_dir = NULL, contour = TRUE){
  group_eemlist <- vector(mode = "list", length = 6)
  class(group_eemlist) <- "eemlist"
  group_eemlist[[1]] <- list_of_eemlists[[1]][[which_eem]]
  group_eemlist[[2]] <- list_of_eemlists[[2]][[which_eem]]
  group_eemlist[[3]] <- list_of_eemlists[[3]][[which_eem]]
  group_eemlist[[4]] <- list_of_eemlists[[4]][[which_eem]]
  group_eemlist[[5]] <- list_of_eemlists[[5]][[which_eem]]
  group_eemlist[[6]] <- list_of_eemlists[[6]][[which_eem]]
  name_2 <- list_of_eemlists[[2]][[which_eem]][["sample"]]
  name_3 <- list_of_eemlists[[3]][[which_eem]][["sample"]]
  name_4 <- list_of_eemlists[[4]][[which_eem]][["sample"]]
  name_5 <- list_of_eemlists[[5]][[which_eem]][["sample"]]
  name_6 <- list_of_eemlists[[6]][[which_eem]][["sample"]]
  group_eemlist[[2]][["sample"]] <- paste0(name_2,"_masked")
  group_eemlist[[3]][["sample"]] <- paste0(name_3,"_interp")
  group_eemlist[[4]][["sample"]] <- paste0(name_4,"_desct")
  group_eemlist[[5]][["sample"]] <- paste0(name_5,"_desct2")
  group_eemlist[[6]][["sample"]] <- paste0(name_6,"_final")
  if(!(is.null(output_dir))){
    png(paste0(output_dir,group_eemlist[[1]][["sample"]],"processing_steps.png"), units="cm", width=18, height=14, res=300)
    print(eem_overview_plot(group_eemlist, contour = contour)) # change here to output contours in .png EEMs
    dev.off()
  } else(
    message("no output directory defined; png not exported")
  )
  print(eem_overview_plot(group_eemlist[1:6], spp = 6, contour = contour))
}
