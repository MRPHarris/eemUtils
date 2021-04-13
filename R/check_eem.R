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

