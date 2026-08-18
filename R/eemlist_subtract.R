#' Subtract an eem from an eemlist
#'
#' @description Take an eem and subtract it from an eemlist. Useful for, e.g., bulk analytical blank subtraction.
#'
#' @param eems_minuend an eemlist compliant with the staRdom/EEM/eemR framework
#' @param eem_subtrahend a single EEM.
#'
#' @export
#'
eemlist_subtract <- function(eems_minuend,
                             eem_subtrahend){
  ## input checks
  if(class(eems_minuend) != "eemlist"){
    stop("object eems_minuend must be of class eemlist")
  }
  if(class(eem_subtrahend) == "eemlist"){
    eem_subtrahend = eem_subtrahend[[1]] # Pull out just the eem.
  }
  ## new eemlist to add processed eems to.
  eemlist_sub <- vector(mode = "list", length = length(eems_minuend))
  class(eemlist_sub) <- "eemlist"
  ## coerce subtrahend to data.frame
  eem_subtrahend_df = as.data.frame(eem_subtrahend, gather = FALSE)
  ## Main subtraction loop
  for(e in seq_along(eems_minuend)){
    ## coerce minuend for this iteration to data frame.
    eem_minuend_it <- eems_minuend[[e]]
    eem_minuend_df_it <- as.data.frame(eem_minuend_it, gather = FALSE)
    ## check that the dimensions are the same.
    if(!isTRUE(nrow(eem_minuend_df_it) == nrow(eem_subtrahend_df) & ncol(eem_minuend_df_it) == ncol(eem_subtrahend_df))){
      stop("EEMs minuend and subtrahend are of differing dimensions.")
    }
    ## FEATURE UPDATE: ADD extend2largest() OPTION HERE
    ## Subtraction
    eem_minuend_df_new <- eem_minuend_df_it - eem_subtrahend_df
    ## Return to eem class object
    eem_minuend_new <- eemdf_to_eem(eemdf = eem_minuend_df_new,
                                    file = eem_minuend_it[['file']],
                                    sample = eem_minuend_it[['sample']],
                                    location = eem_minuend_it[['location']])
    ## Add processed eem from this iteration to the new eemlist
    eemlist_sub[[e]] <- eem_minuend_new
  }
  return(eemlist_sub)
}
