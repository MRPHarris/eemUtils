#' Work in progress function to bin or interpolate series data at a monthly resolution.
#'
#' @description Bin or interpolate data to a monthly scale and output in a dataframe
#'        that is compatible with the KNMI climate explorer user upload field.
#'        Originally developed for use with sub-annually resolved firn core fluorescence
#'        intensity data.
#'
#' @param data The sub-annually resolved series data to be re-scaled.
#' @param old_timescale the original time index for the data.
#' @param new_timescale the new timescale to be re-scaled into. E.g. seq(2010,2000,1/12)
#' @param n_month_groups placeholder - number of months. Leave as 12, as that's how many months there are I guess.
#' @param method 'bin' or 'interpolate' - how will the rescaling be performed? Binning is usually better.
#' @param new_timescale_centered TRUE/FALSE if binning, will the bins be centered?
#' @param replace_NA TRUE/FALSE to replace instances of NA with -999.9, which is KNMI's default missing value handle.
#'
#' @export
#'

knmi_monthly_rescale <- function(data, old_timescale, new_timescale,
                                 n_month_groups = 12, method = 'bin', new_timescale_centered = FALSE, replace_NA = TRUE){
  if(method == "bin"){
    message("Method: Binning")
    xDF <- as.data.frame(data) # make sure inputs can be properly indexed
    Ages <- as.data.frame(old_timescale)
    coreMat = xDF # conform to legacy naming scheme
    if(!isTRUE(new_timescale_centered)){ # Centering bins for better representative accuracy.
      coreBin <- new_timescale+(1/2*(new_timescale[2]-new_timescale[1]))
      bin_lower_edges <- round(coreBin - (1/2*(coreBin[2]-coreBin[1])),2)
      bin_upper_edges <- round(coreBin + (1/2*(coreBin[2]-coreBin[1])),2)
    } else{
      coreBin <- new_timescale
      bin_lower_edges <- round(coreBin - (1/2*(coreBin[2]-coreBin[1])),2)
      bin_upper_edges <- round(coreBin + (1/2*(coreBin[2]-coreBin[1])),2)
    }
    sampleBin = vector("list", length(coreBin))
    # bin widths just in case
    binwidth <- coreBin[2]-coreBin[1]
    binwidth_half <- (1/2*binwidth)
    # binning loop (FOR MONTHLY)
    for(i in 1:length(sampleBin)){
      toBin = which(Ages >= bin_lower_edges[i] & Ages < bin_upper_edges[i]) # currently, this bins between this bin and the next. It is a leading bin, not a centered one.
      sampleBin[[i]] = coreMat[toBin,]
      names(sampleBin)[[i]] <- coreBin[i]
    }
    coreRes = matrix(NA, ncol=ncol(xDF), nrow=length(sampleBin))
    if(ncol(xDF) > 1) { # means of binned data points
      for(i in 1:length(sampleBin)){
        coreRes[i,] = apply(sampleBin[[i]], 2, function(x)mean(x, na.rm=T))
      }
    } else {
      for(i in 1:length(sampleBin)){
        coreRes[i,] = mean(sampleBin[[i]], na.rm=T)
      }
    }
    coreRes_1 <- data.frame(matrix(NA,nrow = nrow(coreRes), ncol = ncol(xDF)+1))
    coreRes_1[,1] <- coreBin
    coreRes_1[,2] <- coreRes
    colnames(coreRes_1) <- c("Age_Bin","values")
    data_scaled <- coreRes_1
    years <- substr(data_scaled[,1], start = 1, stop = 4)
    years_unique <- unique(years)
  }
  if(method == "interp"){ # outputs: new xdata and new ydata
    message("Method: Interpolation with approx(). Be advised that binning is recommended in most use cases")
    if(ncol(as.data.frame(data)) > 1){
      stop("currently no support for >1 columns for interpolation. Instead pass each ydata column to the function separately.")
    }
    # create scale for interpolation.
    if(!isTRUE(new_timescale_centered)){ # Centering bins for better representative accuracy.
      new_scale <- new_timescale+(1/2*(new_timescale[2]-new_timescale[1]))
    } else{
      new_scale <- new_timescale
    }
    # Interpolating here
    data_scaled_list <- approx(x = old_timescale, y = data, xout = new_scale)
    # Tidying
    data_scaled <- data.frame(matrix(NA, nrow = length(data_scaled_list[[1]]),ncol = 2))
    data_scaled[,1] <- round(rev(data_scaled_list[["x"]]),3)
    data_scaled[,2] <- rev(data_scaled_list[["y"]])
    # Now to convert to a useable KNMI format.
    years <- substr(data_scaled[,1], start = 1, stop = 4)
    years_unique <- unique(years)
  }
  # creating frame for transposed data
  knmi_frame <- data.frame(matrix(NA,nrow = length(years_unique), ncol = n_month_groups+1))
  months <- c("","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
  colnames(knmi_frame) <- months # months as column names
  knmi_frame[1:length(years_unique),1] <- years_unique # years into rows
  # loop that will go year-by-year and transpose relevant rows to the above frame.
  it_list <- vector(mode = "list", length = length(years_unique))
  for(i in seq_along(it_list)){
    if(i != 1){
      # one or more indexes for this year in loop
      index_it <- as.numeric(str_which(data_scaled[,1],years_unique[i]))
      # extract rows which match this indexcla
      rows_it <- (data_scaled[index_it,2])
      # transpose these rows into matching row in matrix
      knmi_frame[i,2:(length(rows_it)+1)] <- rows_it
    } else {
      # one or more indexes for this year in loop
      index_it <- as.numeric(str_which(data_scaled[,1],years_unique[i]))
      # extract rows which match this index
      rows_it <- (data_scaled[index_it,2])
      # transpose these rows into matching row in matrix
      knmi_frame[i,(ncol(knmi_frame)-(length(rows_it)-1)):ncol(knmi_frame)] <- rows_it
    }
  }
  if(isTRUE(replace_NA)){
    knmi_frame[is.na(knmi_frame)] <- -999.9 # as per knmi's guidelines on NA handling.
  }
  knmi_frame
}
