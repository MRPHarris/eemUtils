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
  group_eemlist <- vector(mode = "list", length = length(list_of_eemlists))
  class(group_eemlist) <- "eemlist"
  counter = 0
  group_eemlist <- lapply(list_of_eemlists, function(eem){
    counter = counter + 1
    list_of_eemlists[[counter]][[which_eem]]
  })


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
