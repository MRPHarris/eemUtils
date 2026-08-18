#' Sum eems within an eemlist together.
#'
#' @description Simple addition of the EEM objects within a supplied eemlist.
#'
#' @param eemlist A list of EEMs in a format compliant with eemR/staRdom.
#'
#' @export
#'
eemlist_sum <- function(eemlist){
  # Input check
  if(class(eemlist) != 'eemlist'){
    stop("Please pass the function an object of class 'eemlist'")
  }
  # Construct data frames
  eemlist_dataframes <- lapply(eemlist,function(e){
    e_n <- as.data.frame(e, gather = FALSE)
    e_n
  })
  # Perform addition
  adds <- Reduce("+", eemlist_dataframes)
  adds <- eemdf_to_eem(eemdf = adds,
                       file = NULL,
                       sample = 'summed_eemlist',
                       location = NULL,
                       gathered = FALSE)
  adds
}
