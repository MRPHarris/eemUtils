#' Change wavelength range of eems with in an eemlist.
#'
#' @description A tweaked/re-built version of staRdom::eem_range(), which was not
#'      working for some datasets. Use this function to change the wavelength range
#'      of a group of EEMs.
#'
#' @param eemlist A list of EEMs compliant with the eemR/staRdom framework.
#' @param ex_range Numeric excitation min:max of the new range to constrict the EEMs to.
#' @param em_range Numeric excitation min:max of the new range to constrict the EEMs to.
#'
#' @export
#'

eem_range_mod <- function(eemlist, ex_range, em_range){
  if(class(eemlist) != "eemlist"){
    warning("please pass the function an object of class eemlist")
  }
  eemlist_name <- deparse(substitute(eemlist))
  # establish limits by indexing first EEM.
  # excitation
  min_ex <- as.numeric(min(ex_range)) #get min/max
  max_ex <- as.numeric(max(ex_range))
  ex_nm <- eemlist[[1]]$ex    # extract exitation wavelengths
  ex_nm <- ex_nm[length(ex_nm):1]   # flip to allow indexing
  dt_ex <- data.table(ex_nm, val = ex_nm)   # data.table
  setattr(dt_ex,"sorted", "ex_nm")   # set data.table attributes so it can be indexed
  min_ex_index_r <- as.numeric(dt_ex[J(max_ex), roll = "nearest", which = TRUE]) # get min index
  min_ex_index <- length(eemlist[[1]]$ex) - (min_ex_index_r - 1)   # flip, as the true wavelengths are the other way round
  max_ex_index_r <- as.numeric(dt_ex[J(min_ex), roll = "nearest", which = TRUE]) # same as above, for max
  max_ex_index <- length(eemlist[[1]]$ex) - (max_ex_index_r - 1)
  # emission
  min_em <- as.numeric(min(em_range))
  max_em <- as.numeric(max(em_range))
  em_nm <- eemlist[[1]]$em   # acquiring nearest value in table.
  dt_em <- data.table(em_nm, val = em_nm)
  setattr(dt_em, "sorted", "em_nm")
  min_em_index <- as.numeric(dt_em[J(min_em), roll = "nearest", which = TRUE]) # indem of said value
  max_em_index <- as.numeric(dt_em[J(max_em), roll = "nearest", which = TRUE])
  is_between <- function(x, a, b) {
    x >= a & x <= b
  }
  new_eemlist <- eemlist
  for (i in seq_along(eemlist)){
    # excitation range
    if (!missing(ex_range)){
      stopifnot(is.numeric(ex_range), all(ex_range >= 0))
      index_ex <- which(is_between(eemlist[[i]]$ex, eemlist[[i]]$ex[max_ex_index],
                                   eemlist[[i]]$ex[min_ex_index]))
      if (length(index_ex != 0)) {
        eemlist[[i]]$ex <- eemlist[[i]]$ex[index_ex]
        eemlist[[i]]$x <- eemlist[[i]]$x[, index_ex]
      }
    }
    # emission range
    if (!missing(em_range)){
      stopifnot(is.numeric(em_range), all(em_range >= 0))
      index_em <- which(is_between(eemlist[[i]]$em, eemlist[[i]]$em[min_em_index],
                                   eemlist[[i]]$em[max_em_index]))
      if (length(index_em != 0)) {
        eemlist[[i]]$em <- eemlist[[i]]$em[index_em]
        eemlist[[i]]$x <- eemlist[[i]]$x[index_em, ]
      }
    }
    new_eemlist[[i]] <- eemlist[[i]]
  }
  new_eemlist
}
