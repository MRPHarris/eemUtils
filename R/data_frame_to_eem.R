#' Attempt to convert a raw EEM data dataframe to an eemR-compliant object of class 'eem'.
#'
#' @description Takes a data frame (e.g. raw EEM data) and attempts to coerce it into an eem object compliant with the eemR/staRdom framework.
#'
#' @param x The target dataframe.
#' @param sample Sample/EEM name.
#' @param location Optional; where is the EEM from
#' @param filename Full file name and path of the EEM.
#' @param ex Which side is the excitation axis on - columns, or rows?
#'
#' @export
#'

data_frame_to_eem <- function(x, sample = NULL, filename = NULL, location = NULL, ex = 'cols'){
  if(!(ex == 'cols' || ex == 'rows')){
    stop("ex must be either 'cols' or 'rows'")
  }
  eem <- vector(mode = 'list', length = 6)
  names(eem) <- c("file","sample","x","ex","em","location")
  # 6 parts to the list.
  eem[["file"]] <- filename
  eem[["location"]] <- location
  eem[["sample"]] <- sample
  eem[["x"]] <- as.matrix(x)
  if(ex == 'cols'){
    eem[["ex"]] <- as.numeric(colnames(x))
    eem[["em"]] <- as.numeric(rownames(x))
  } else{
    eem[["ex"]] <- as.numeric(rownames(x))
    eem[["em"]] <- as.numeric(colnames(x))
  }
  e <- function(y){structure(y,class = 'eem')}
  eem <- e(eem)
}
