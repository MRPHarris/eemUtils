#' Use data.table to find the nearest value in a vector.
#'
#' @description Using binary search provided by the data table package, find the nearest value.
#'
#' @param data a vector of numeric values to search within
#' @param value to drive binary search. i.e. output will be closest value to this input.
#'
#' @export
#'
binary_search_nearest <- function(data, value){
  if(!is.numeric(data)){
    data <- as.numeric(data)
  }
  if(!is.numeric(value)){
    value <- is.numeric(value)
  }
  data_tab <- data.table(data, val = data)
  setattr(data_tab, "sorted","data")
  setkey(data_tab, data)
  nearest <- as.numeric(data_tab[J(value), roll = "nearest"][1,2])
  nearest
}
