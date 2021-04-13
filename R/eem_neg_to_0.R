#' Set all negative values within a group of EEMs to 0.
#'
#' @description Negative fluorescence is not possible, and typically indicates
#'        bad processing, noise, or artefacts. This function will set all cases
#'        of negative fluorescence to 0.
#'
#' @param eemlist A list of EEMs in a format compliant with eemR/staRdom.
#' @param outputfolder Path to the folder where the new eemlist will be sent. Or NULL, if no export is desired.
#'
#' @export
#'

# This function removes sets all negative values to 0.
eem_neg_to_0 <- function(eemlist, outputfolder = NULL){
  EEMs_NoNeg <- eemlist
  for(i in seq_along(EEMs_NoNeg)){                                         # main for loop
    eem_ungathered <- as.data.frame(EEMs_NoNeg[[i]], gather = FALSE)       # extract EEM, don't gather
    eem_ungathered[,][eem_ungathered[,] <0] <- 0                        # set all values less than 0 in EEM to 0
    if(!is.null(outputfolder)){
      write.csv(eem_ungathered, file = paste0(outputfolder,EEMs_NoNeg[[i]][["sample"]],"_noneg.csv"), row.names = TRUE) # Export EEM with iterative naming scheme.
    }
  }
  EEMs_NoNeg <- EEMs_NoNeg
}
