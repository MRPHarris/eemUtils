#' Set areas of one or more EEMs to NA.
#'
#' @description A tweaked/re-built version of staRdom::set_NA(). Sets areas of
#'        eems to NA. Option to interpolate. Changes from set_NA() relate to how
#'        the range is identified. I recommend interpolating seperately.
#'
#' @param eemlist A list of EEMs in a format compliant with eemR/staRdom.
#' @param sample optional, names or indices of samples to process
#' @param em optional, emission wavelengths to set NA
#' @param ex optional, excitation wavelengths to set NA
#' @param interpolate FALSE, 1 or 2, interpolate NAs or not, 2 different methods.
#'
#' @export
#'

eem_setNA_mod <- function(eem_list, sample = NULL, em = NULL, ex = NULL, interpolate = FALSE){
  if (is.null(sample)) {
    sample <- eem_names(eem_list)
  }
  if (is.numeric(sample)) {
    sample <- eem_names(eem_list)[sample]
  }
  eem_list <- lapply(eem_list, function(eem) {
    if (eem$sample %in% sample) {
      if (is.null(ex)){
        ex2 <- 1:ncol(eem$x)
      } else (ex2 <- which(eem$ex %in% ex))
      if (is.null(em)){
        em2 <- 1:nrow(eem$x)
      }
      if(length(em) > 1){
        min_em <- as.numeric(min(em))
        max_em <- as.numeric(max(em))
        em_nm <- eem$em   # acquiring nearest value in table.
        dt_em <- data.table(em_nm, val = em_nm)
        setattr(dt_em, "sorted", "em_nm")
        min_em_val <- as.numeric(dt_em[J(min_em), roll = "nearest"][1,2]) # Identify value closest to that specified
        min_em_index <- as.numeric(dt_em[J(min_em), roll = "nearest", which = TRUE]) # index of said value
        max_em_val <- as.numeric(dt_em[J(max_em), roll = "nearest"][1,2])
        max_em_index <- as.numeric(dt_em[J(max_em), roll = "nearest", which = TRUE])
        em <- eem$em[min_em_index:max_em_index]
        em2 <- which(eem$em %in% em)}  # new em comprised of actual values
      if(length(em) == 1){
        em_nm <- eem$em
        dt_em <- data.table(em_nm, val = em_nm)
        setattr(dt_em, "sorted", "em_nm")
        em_val <- as.numeric(dt_em[J(em), roll = "nearest"][1,2])
        em_index <- as.numeric(dt_em[J(em), roll = "nearest", which = TRUE])
        em2 <- which(eem$em %in% em)}
      eem$x[em2, ex2] <- NA
    }
    eem
  }) %>% `class<-`("eemlist")
  if (interpolate != FALSE) {
    eem_list[which(eem_names(eem_list) %in% sample)] <- eem_interp(eem_list[which(eem_names(eem_list) %in%
                                                                                    sample)], type = interpolate) %>% `class<-`("eemlist")
  }
  eem_list
}
