#' A work in progress function to create mirrored offsets of series data.
#'
#' @description # Extracts xdata (an index e.g. depth) and ydata and creates
#'        mirrored offsets, such that the lines between bars in a ggplot2 geom_bar
#'        plot will be obscured if set to white.At the moment, columns must match
#'        the default. yoffset can be changed to suit ydata range. This function
#'        is very much a WIP.
#'
#' @param Input_Data_df The input data frame. Should comprise of at least 1 index (x axis) and the variable (y axis).
#' @param xdata_column column for the series variable (x-axis)
#' @param data_column column containing the independent/y-axis variable.
#' @param yoffset adjustable offset. Use in conjunction with plotting code.
#'
#' @export
#'

barify <- function(Input_Data_df, xdata_column = 1, data_column = 2, yoffset = 0.5){
  colname_xdata <- colnames(Input_Data_df[c(xdata_column)])
  colname_data <- colnames(Input_Data_df[c(data_column)])
  data <- as.numeric(as.matrix(Input_Data_df[data_column]))
  xdata <- as.numeric(as.matrix(Input_Data_df[xdata_column]))
  Input_Data <- cbind(as.data.frame(xdata),as.data.frame(data))
  DataDiff <- Input_Data %>%
    mutate(diff_data = data - lag(data, default = first(data)))
  DataDiff <- DataDiff %>%
    mutate(diff_xdata = xdata - lag(xdata, default = first(xdata)))
  DataDiff <- DataDiff
  DataDiff2 <- DataDiff
  shift <- function(x, n){
    c(x[-(seq(n))], rep(0, n))
  }
  DataDiff$diff_data <- shift(DataDiff$diff_data,1)
  DataDiff <- within(DataDiff, data[diff_data<0] <- NA)
  DataDiff <- within(DataDiff, xdata <- xdata + 0.5*(abs(diff_xdata)))
  DataDiff <- within(DataDiff, data[!is.na(diff_data)] <- (data[!is.na(diff_data)] - yoffset))
  DataDiff2 <- within(DataDiff2, data[!diff_data <= 0] <- NA)
  DataDiff2 <- within(DataDiff2, xdata <- xdata - 0.5*(abs(diff_xdata)))
  DataDiff2 <- within(DataDiff2, data[!is.na(diff_data)] <- (data[!is.na(diff_data)] - yoffset))
  DataDiff2 <- within(DataDiff2, xdata[1] <- (xdata[1] - 0.5*abs(diff_xdata[2])))
  colnames(DataDiff)[data_column] <- colname_data
  colnames(DataDiff)[xdata_column] <- colname_xdata
  colnames(DataDiff2)[data_column] <- colname_data
  colnames(DataDiff2)[xdata_column] <- colname_xdata
  DataDiff <<- DataDiff
  DataDiff2 <<- DataDiff2
  assign(paste0(deparse(substitute(Input_Data_df)),"_OffsetA"), DataDiff, envir = parent.frame())
  assign(paste0(deparse(substitute(Input_Data_df)),"_OffsetB"), DataDiff2, envir = parent.frame())
  rm(DataDiff,DataDiff2,envir = parent.frame())
}
