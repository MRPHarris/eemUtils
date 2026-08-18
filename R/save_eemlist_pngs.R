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
    class(Target_eemlist) <- "eemlist"
    Target_eemlist[[1]] <- eemlist[[i]]
    png(paste0(output_dir,EEM_Name,".png"), units="cm", width=18, height=14, res=300)
    print(ggeem2(Target_eemlist, contour = contour)) # change here to output contours in .png EEMs
    dev.off()
  }
}
