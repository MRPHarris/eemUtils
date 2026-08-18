#' Average a set of EEMs.
#'
#' @description Average the EEMs within an eemlist.
#'
#' @param eemlist A list of EEMs, compliant with the eemR/staRdom framework.'
#' @param verbose TRUE/FALSEto show message at completion of operation.
#'
#' @export
#'
eemlist_average <- function(eemlist, verbose = TRUE){
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
    if(isTRUE(verbose)){
      message(paste0(length(eemlist)," EEMs passed to average_eems() for averaging."))
    }
    new_eemlist[[1]] <- averaged_eem
    return(new_eemlist)
  }
}
