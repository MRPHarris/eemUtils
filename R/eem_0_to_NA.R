#' Set all instances of '0'/zero to NA within a list of EEMs.
#'
#' @description Go through a list of eemR/staRdom-compliant EEMs, and replace all
#'       instances of 0 with NA. Used in cases where e.g. an Aqualog Rayleigh-masking
#'       step has masked those rows with 0, preventing proper interpolation and
#'       PARAFAC modeling.
#'
#' @param eemlist A list of EEMs in a format compliant with eemR/staRdom.
#' @param outputfolder Path to the folder where the new eemlist will be sent. Or NULL, if no export is desired.
#'
#' @export
#'

eem_0_to_NA <- function(eemlist, outputfolder = NULL){
  EEMs_DeNeg <- eemlist# load required packages
  for(i in seq_along(EEMs_DeNeg)){                                         # main for loop
    eem_ungathered <- as.data.frame(EEMs_DeNeg[[i]], gather = FALSE)       # extract EEM, don't gather
    eem_ungathered[,][eem_ungathered[,] == 0] <- NA                        # set all values of 0 in EEM to NA
    if(!is.null(outputfolder)){
      write.csv(eem_ungathered, file = paste0(outputfolder,EEMs_DeNeg[[i]][["sample"]],"_masked.csv"), row.names = TRUE) # Export EEM with iterative naming scheme.
    }
  }
  EEMs_DeNeg <- EEMs_DeNeg  # re-imports eems from the folder they were exported to
}
