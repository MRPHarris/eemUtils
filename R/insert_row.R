#' Insert a single row seamlessly into a data frame.
#'
#' @description Insert a row into a data frame whilst shifting the other rows to
#'      make room.
#'
#' @param dataframe The target data frame.
#' @param newrow_index numeric - row position will the new row occupy?
#' @param newrow_contents The new row to be inserted into the dataframe.
#'
#' @export
#'

# Insert a row into a dataset, e.g. for missing values.
insert_row <- function(dataframe, newrow_index, newrow_contents){
  nrow_init <- as.numeric(nrow(dataframe))
  ncol_init <- ncol(dataframe)-1
  insert_frame <- data.frame(matrix(nrow = nrow_init+1, ncol = ncol_init), data = NA)
  colnames(insert_frame) <- colnames(dataframe)
  indexA <- newrow_index-1
  indexB <- newrow_index+1
  insert_frame[1:indexA,] <- dataframe[1:indexA,]
  insert_frame[indexB:nrow(insert_frame),] <- dataframe[newrow_index:nrow(dataframe),]
  insert_frame[newrow_index,] <- newrow_contents
  insert_frame
}

