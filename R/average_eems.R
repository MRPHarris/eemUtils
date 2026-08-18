#' Average the eems within an EEMlist,
#'
#' @description Average all the pixels within an eemlist, producing a per-ex/em pixel average EEM object
#'
#' @param eemlist an eemlist compliant with the staRdom/EEM/eemR framework
#'
#' @export
#'
average_eems <- function(eemlist){
  if(length(eemlist) == 1){
    message("1 EEM passed to average_eems() for averaging. Returning unchanged.")
    new_eemlist <- eemlist
    return(new_eemlist)
  } else if(length(eemlist) > 1){
    ungathered_list <- vector(mode = "list",length = length(eemlist))
    for(i in seq_along(eemlist)){
      eem_it <- eemlist[[i]]
      file_it <- eem_it[['file']]
      sample_it <- eem_it[['sample']]
      location_it <- eem_it[['location']]
      eem_ungathered <- as.data.frame(eemlist[[i]], gather = FALSE)       # extract EEM, don't gather
      ungathered_list[[i]] <- eem_ungathered
    }
    # Average it
    averaged <- Reduce("+",ungathered_list)/length(ungathered_list)
    # Now back to eem
    averaged_eem <- eemdf_to_eem(averaged,
                                 file = "",
                                 sample = "averaged_eem",
                                 location = "")
    new_eemlist <- vector(mode = "list",length = 1)
    class(new_eemlist) <- "eemlist"
    message(paste0(length(eemlist)," EEMs passed to average_eems() for averaging."))
    new_eemlist[[1]] <- averaged_eem
    return(new_eemlist)
  }

}
