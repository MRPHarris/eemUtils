#' Save a list of EEMs as individual .csv files.
#'
#' @description # Save the EEMs comprising an eemlist to a set of csv files. If no
#'      output_dir is set, the csvs are written straight to the working directory.
#'
#' @param eemlist A list of EEMs in a format compliant with eemR/staRdom.
#' @param outputfolder Either the path of the folder to write the .csvs to, or NULL for the working directory.
#' @param append_name optional; a character vector that will be appended to each .csv file.
#'
#' @export
#'

save_eemlist_csvs <- function(eemlist,outputfolder = NULL, append_name = NULL){
  if(is.null(outputfolder)){
    output_dir <- getwd()
    message("No outputfolder specified; exporting to working directory")
  }
  if(!is.null(append_name)){
    append_name = paste0("_",append_name)
  }
  for (i in seq_along(eemlist)){
    eem_ungathered <- as.data.frame(eemlist[[i]], gather = FALSE)
    write.csv(eem_ungathered, file = paste0(outputfolder,eemlist[[i]][["sample"]],append_name,".csv"))
    message(i,"/",length(eemlist)," | ",format(round((i/length(eemlist))*100,1), nsmall = 1),"% complete")
  }
  message("Done!")
}
