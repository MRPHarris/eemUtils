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
#' @param annual_average_startmonth NULL or numeric between 1 and 12. If set to a number between 1 and 12, an annual average will be calculated starting at that month.
#' @param average_type one of either "annual_averages" or "impersonate_monthly".. The former will return the average values in data frame, and the latter will impersonate a format that KNMI will accept as a monthly resolution data to allow comparison of the annual averages against monthly-resolved data.
#' @param average_alignment NULL or one of either "current" or "previous". NULL means that Jan to June will (months 1 to 6) will commence the average in the current year, and July to Dec in the preceding year.
#' @param align_monthly_impersonation TRUE/FALSE. If returning an impersonate_monthly data frame (see previous param), should the annual values start in January for a given year, or commence at the annual_average_startmonth?
#'
#' @export
#'
knmi_monthly_rescale <- function(data,
                                 old_timescale,
                                 new_timescale,
                                 n_month_groups = 12,
                                 method = 'bin',
                                 new_timescale_centered = FALSE,
                                 replace_NA = TRUE,
                                 annual_average_startmonth = NULL,
                                 average_type = "impersonate_monthly",
                                 average_alignment = NULL,
                                 align_monthly_impersonation = FALSE){
  if (method == "bin") {
    message("Method: Binning")
    xDF <- as.data.frame(data)
    Ages <- as.data.frame(old_timescale)
    coreMat = xDF
    if (!isTRUE(new_timescale_centered)) {
      coreBin <- new_timescale + (1/2 * (new_timescale[2] -
                                           new_timescale[1]))
      bin_lower_edges <- round(coreBin - (1/2 * (coreBin[2] -
                                                   coreBin[1])), 2)
      bin_upper_edges <- round(coreBin + (1/2 * (coreBin[2] -
                                                   coreBin[1])), 2)
    } else {
      coreBin <- new_timescale
      bin_lower_edges <- round(coreBin - (1/2 * (coreBin[2] -
                                                   coreBin[1])), 2)
      bin_upper_edges <- round(coreBin + (1/2 * (coreBin[2] -
                                                   coreBin[1])), 2)
    }
    sampleBin = vector("list", length(coreBin))
    binwidth <- coreBin[2] - coreBin[1]
    binwidth_half <- (1/2 * binwidth)
    for (i in 1:length(sampleBin)) {
      toBin = which(Ages >= bin_lower_edges[i] & Ages <
                      bin_upper_edges[i])
      sampleBin[[i]] = coreMat[toBin, ]
      names(sampleBin)[[i]] <- coreBin[i]
    }
    coreRes = matrix(NA, ncol = ncol(xDF), nrow = length(sampleBin))
    if (ncol(xDF) > 1) {
      for (i in 1:length(sampleBin)) {
        coreRes[i, ] = apply(sampleBin[[i]], 2, function(x) mean(x,
                                                                 na.rm = T))
      }
    } else {
      for (i in 1:length(sampleBin)) {
        coreRes[i, ] = mean(sampleBin[[i]], na.rm = T)
      }
    }
    coreRes_1 <- data.frame(matrix(NA, nrow = nrow(coreRes),
                                   ncol = ncol(xDF) + 1))
    coreRes_1[, 1] <- coreBin
    coreRes_1[, 2] <- coreRes
    colnames(coreRes_1) <- c("Age_Bin", "values")
    data_scaled <- coreRes_1
    years <- substr(data_scaled[, 1], start = 1, stop = 4)
    years_unique <- unique(years)
  }
  if (method == "interp") {
    message("Method: Interpolation with approx(). Be advised that binning is recommended in most use cases")
    if (ncol(as.data.frame(data)) > 1) {
      stop("currently no support for >1 columns for interpolation. Instead pass each ydata column to the function separately.")
    }
    if (!isTRUE(new_timescale_centered)) {
      new_scale <- new_timescale + (1/2 * (new_timescale[2] -
                                             new_timescale[1]))
    }
    else {
      new_scale <- new_timescale
    }
    data_scaled_list <- approx(x = old_timescale, y = data,
                               xout = new_scale)
    data_scaled <- data.frame(matrix(NA, nrow = length(data_scaled_list[[1]]),
                                     ncol = 2))
    data_scaled[, 1] <- round(rev(data_scaled_list[["x"]]),
                              3)
    data_scaled[, 2] <- rev(data_scaled_list[["y"]])
    years <- substr(data_scaled[, 1], start = 1, stop = 4)
    years_unique <- unique(years)
  }
  knmi_frame <- data.frame(matrix(NA, nrow = length(years_unique),
                                  ncol = n_month_groups + 1))
  months <- c("",seq(1,n_month_groups,1))
  colnames(knmi_frame) <- months
  knmi_frame[1:length(years_unique), 1] <- years_unique
  it_list <- vector(mode = "list", length = length(years_unique))
  for (i in seq_along(it_list)) {
    if (i != 1) {
      index_it <- as.numeric(str_which(data_scaled[, 1],
                                       years_unique[i]))
      rows_it <- (data_scaled[index_it, 2])
      knmi_frame[i, 2:(length(rows_it) + 1)] <- rows_it
    }
    else {
      index_it <- as.numeric(str_which(data_scaled[, 1],
                                       years_unique[i]))
      rows_it <- (data_scaled[index_it, 2])
      knmi_frame[i, (ncol(knmi_frame) - (length(rows_it) -
                                           1)):ncol(knmi_frame)] <- rows_it
    }
  }
  # if (isTRUE(replace_NA)) {
  #   knmi_frame[is.na(knmi_frame)] <- -999.9
  # }
  if(!is.null(annual_average_startmonth)){
    # Computing annual averages
    knmi_frame_rev <- knmi_frame %>%
      'colnames<-'(c("year",colnames(knmi_frame)[2:ncol(knmi_frame)])) %>%
      column_to_rownames(var = "year")
    if(average_type == "impersonate_monthly"){
      knmi_frame_avs <- data.frame(matrix(NA,nrow = nrow(knmi_frame_rev), ncol = ncol(knmi_frame_rev))) %>%
        'colnames<-'(colnames(knmi_frame_rev)) %>%
        'rownames<-'(rownames(knmi_frame_rev))
    }
    # For a given starting month,
    ann_avs <- data.frame(matrix(NA,nrow = nrow(knmi_frame_rev), ncol = 3)) %>%
      'colnames<-'(c('year','value','count'))
    ann_avs$year <- as.numeric(rownames(knmi_frame_rev))
    # making a loop because I can't wrap my head around this at the moment
    # how many years?
    n_years <- nrow(knmi_frame_rev)
    it_list <- vector('list', length = n_years)
    # monthly alignment added. Default behaviour by setting 'average_alignment' to NULL.
    if(is.null(average_alignment)){
      if(annual_average_startmonth <= 6){
        type = 'current'
        message("Annual average will start in the current calendar year.")
      } else if(annual_average_startmonth >= 7){
        type = 'previous'
        message("Annual average will start in the preceding calendar year.")
      }
    } else if(!is.null(average_alignment)){
      if(average_alignment != "current" && average_alignment != "previous"){
        stop("average alignment must be one of NULL, 'current' or 'previous'.")
      }
      type = average_alignment
    }
    # if(annual_average_startmonth > 12){
    #   stop("There are only 12 months in a year. annual_average_startmonth must be between or equal to 1 and 12.")
    # } else if(annual_average_startmonth <= 6){
    #   type = 'current'
    #   message("Annual average will start in the current calendar year.")
    # } else if(annual_average_startmonth >= 7){
    #   type = 'previous'
    #   message("Annual average will start in the preceding calendar year.")
    # }
    smonth_column <- as.numeric(which(as.numeric(colnames(knmi_frame_rev)) == annual_average_startmonth))
    for(y in seq_along(it_list)){
      # Determine position of starting month in previous year
      year_it <- ann_avs$year[y]
      # Where is this year in the knmi frame?
      year_row <- which(ann_avs$year == year_it)
      # get average value for this year
      if(type == "previous"){
        prev_row_vals <- as.numeric(knmi_frame_rev[year_row-1,])[seq(smonth_column,ncol(knmi_frame_rev),1)]
        current_row_vals <- as.numeric(knmi_frame_rev[year_row,])[seq(1,smonth_column-1,1)]
        vals <- c(prev_row_vals,current_row_vals)
        ann_avs$value[year_row] = mean(vals, na.rm = T)
        ann_avs$count[year_row] = sum(!is.na(vals))
        if((average_type == "impersonate_monthly") && (isTRUE(align_monthly_impersonation))){
          if(y != 1){
            knmi_frame_avs[year_row-1,][seq(smonth_column,ncol(knmi_frame_rev),1)] <- mean(vals, na.rm = T)
          }
          knmi_frame_avs[year_row,][seq(1,smonth_column-1,1)] <- mean(vals, na.rm = T)
        }
      } else if(type == 'current'){
        current_row_vals <- as.numeric(knmi_frame_rev[year_row,])[seq(smonth_column,ncol(knmi_frame_rev),1)]
        if(smonth_column == 1){
          following_row_vals <- c()
        } else {
          following_row_vals <- as.numeric(knmi_frame_rev[year_row + 1,])[seq(1,smonth_column-1,1)]
        }
        vals <- c(current_row_vals,following_row_vals)
        ann_avs$value[year_row] = mean(vals, na.rm = T)
        ann_avs$count[year_row] = sum(!is.na(vals))
        if((average_type == "impersonate_monthly") && (isTRUE(align_monthly_impersonation))){
          knmi_frame_avs[year_row,][seq(smonth_column,ncol(knmi_frame_rev),1)] <- mean(vals, na.rm = T)
          knmi_frame_avs[year_row + 1,][seq(1,smonth_column-1,1)] <- mean(vals, na.rm = T)
        }
      }
    }
    if((average_type == "impersonate_monthly") && (!isTRUE(align_monthly_impersonation))){
      knmi_frame_avs[,1:12] <- ann_avs$value
    }
    if(average_type == "impersonate_monthly"){
      knmi_frame_avs <- knmi_frame_avs %>%
        rownames_to_column(" ")
      if(n_month_groups == 12){
        months <- c("", "Jan", "Feb", "Mar",
                    "Apr", "May", "Jun", "Jul", "Aug",
                    "Sep", "Oct", "Nov", "Dec")
      }
      if (isTRUE(replace_NA)) {
        knmi_frame_avs[is.na(knmi_frame_avs)] <- -999.9
      }
      return(knmi_frame_avs)
    } else if(average_type == "annual_averages"){
      return(ann_avs)
    }
  } else {
    # Don't undertake any annual average shenanigans
    if(n_month_groups == 12){
      months <- c("", "Jan", "Feb", "Mar",
                  "Apr", "May", "Jun", "Jul", "Aug",
                  "Sep", "Oct", "Nov", "Dec")
    }
    colnames(knmi_frame) <- months
    if (isTRUE(replace_NA)) {
      knmi_frame[is.na(knmi_frame)] <- -999.9
    }
    return(knmi_frame)
  }
}
