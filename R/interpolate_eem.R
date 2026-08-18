#' Add data points between existing values across a specified EEM axis.
#'
#' @description Interpolate between existing eem data points, adding values across a given axis.
#'     Intended for data to be plotted, as an alternative to staRdom's eem_smooth() smoothing function.
#'     Support only for excitation axis at this stage. Not to be confused with staRdom's eem_interp.
#'
#' @param eem An eem or eemlist object compliant with the staRdom/eemR framework.
#' @param n_pp numeric; how many data points to add between existing ones?
#' @param Direction "ex" only.
#'
#' @export
#'
interpolate_eem <- function(eem, n_pp = 2, direction = "ex", verbose = FALSE){
  #eem = DSS0506_MBT_eems_avg[[1]]
  #n_pp = 2
  #direction = "ex"
  if(class(eem) == 'eemlist'){
    eem2 <- lapply(eem, function(e){interpolate_eem(e)}) %>% 'class<-'(c('eemlist'))
    return(eem2)
  } else if(class(eem) == 'eem'){
    eem_df_ug <- as.data.frame(eem, gather = FALSE)
    #eem_df_ug_tofill <- data.frame(matrix(NA, nrow = nrow(eem_df_ug), ncol = ncol(eem_df_ug)))
    #rownames(eem_df_ug_tofill) <- rownames(eem_df_ug)
    #colnames(eem_df_ug_tofill) <- colnames(eem_df_ug)
    if(direction == "ex"){
      # create new EEM frame to fill.
      eem_df_ug_tofill <- data.frame(matrix(NA, nrow = nrow(eem_df_ug), ncol = (ncol(eem_df_ug)*n_pp)-ceiling(n_pp/2)))
      rownames(eem_df_ug_tofill) <- rownames(eem_df_ug)
      min_ex <- as.numeric(min(colnames(eem_df_ug)))
      max_ex <- as.numeric(max(colnames(eem_df_ug)))
      # Iterate along emission axis, deriving excitation scans.
      it_list <- vector(mode = "list", length = nrow(eem_df_ug))
      for(i in seq_along(it_list)){
        ## Derive ex scan for this iteration
        # ex scan
        ex_scan_it <- as.numeric(eem_df_ug[i,])
        ex_scan_it_exvals <- colnames(eem_df_ug)
        ex_scan_it_df <- as.data.frame(matrix(NA, ncol = 2, nrow = length(ex_scan_it)))
        colnames(ex_scan_it_df) <- c("x","y")
        #ex_scan_it_df$x <- seq(1,length(ex_scan_it),1)
        ex_scan_it_df$x <- ex_scan_it_exvals
        ex_scan_it_df$y <- ex_scan_it
        #plot(ex_scan_it_df, type = 'p')
        ## Interpolate this ex scan.
        ex_width <- as.numeric(ex_scan_it_df$x[2]) - as.numeric(ex_scan_it_df$x[1])
        ex_scan_it_extrap <- as.data.frame(approx(x = ex_scan_it_df$x,
                                                  y = ex_scan_it_df$y,
                                                  n = (length(ex_scan_it)*n_pp)-1,
                                                  xout = seq(min(ex_scan_it_df$x),max(ex_scan_it_df$x), ex_width/n_pp)))
        #plot(ex_scan_it_extrap, type = 'p')
        ## Where should the NAs be?
        ## Direct matching
        NA_rows <- as.numeric(ex_scan_it_df$x[which(is.na(ex_scan_it_df$y))])
        # Are there NA rows?
        if(length(NA_rows) > 0){
          # are there gaps?
          gaps <- (NA_rows - lag(NA_rows, default = 0))[-1]
          gap <- which(gaps > ex_width)+1
          if(length(gap) > 0){
            # There are multiple NA areas. 2 max; no support for non-RM missing areas.
            list_it <- vector(mode = "list", length = length(gap)+1)
            # Iterate through each NA section.
            for(g in seq_along(list_it)){
              if(g == 1){
                NA_rows_a <- NA_rows[1:gap-1]
                # first iteration
                NA_rows_start <- as.numeric(which(ex_scan_it_extrap$x == as.numeric(min(NA_rows_a))))
                NA_rows_end <- as.numeric(which(ex_scan_it_extrap$x == as.numeric(max(NA_rows_a))))
                # Need to expand this by half the n_pp value to ensure the rayleigh NA remain as before.
                if(NA_rows_start == 1){
                  NA_rows_end <- NA_rows_end + ceiling(n_pp/2)
                } else {
                  NA_rows_start <- NA_rows_start - ceiling(n_pp/2)
                  NA_rows_end <- NA_rows_end + ceiling(n_pp/2)
                }
                # generate these row indices
                rows <- seq(NA_rows_start,NA_rows_end,1)
                # set these rows to NA
                ex_scan_it_extrap$y[rows] <- NA
              } else {
                NA_rows_a <- NA_rows[gap:length(NA_rows)]
                # first iteration
                NA_rows_start <- as.numeric(which(ex_scan_it_extrap$x == as.numeric(min(NA_rows_a))))
                NA_rows_end <- as.numeric(which(ex_scan_it_extrap$x == as.numeric(max(NA_rows_a))))
                # Need to expand this by half the n_pp value to ensure the rayleigh NA remain as before.
                if(NA_rows_start == 1){
                  NA_rows_end <- NA_rows_end + ceiling(n_pp/2)
                } else {
                  # what if the last NA row is the end of the ex scan?
                  if(ex_scan_it_extrap$x[NA_rows_end] == max_ex){
                    NA_rows_start <- NA_rows_start - ceiling(n_pp/2)
                  } else {
                    NA_rows_start <- NA_rows_start - ceiling(n_pp/2)
                    NA_rows_end <- NA_rows_end + ceiling(n_pp/2)
                  }
                }
              }
              # generate these row indices
              rows <- seq(NA_rows_start,NA_rows_end,1)
              # set these rows to NA
              ex_scan_it_extrap$y[rows] <- NA
            }
          } else {
            ## Where are the NA's going to end up in the extrapolated dataset?
            NA_rows_start <- as.numeric(which(ex_scan_it_extrap$x == as.numeric(min(NA_rows))))
            NA_rows_end <- as.numeric(which(ex_scan_it_extrap$x == as.numeric(max(NA_rows))))
            # Need to expand this by half the n_pp value to ensure the rayleigh NA remain as before.
            if(NA_rows_start == 1){
              NA_rows_end <- NA_rows_end + ceiling(n_pp/2)
            } else {
              # what if the last NA row is the end of the ex scan?
              if(ex_scan_it_extrap$x[NA_rows_end] == max_ex){
                NA_rows_start <- NA_rows_start - ceiling(n_pp/2)
              } else {
                NA_rows_start <- NA_rows_start - ceiling(n_pp/2)
                NA_rows_end <- NA_rows_end + ceiling(n_pp/2)
              }
            }
            # generate these row indices
            rows <- seq(NA_rows_start,NA_rows_end,1)
            # set these rows to NA
            ex_scan_it_extrap$y[rows] <- NA
          }
        }
        # Add extrapolated scan to the frame
        eem_df_ug_tofill[i,] <- ex_scan_it_extrap$y
        if(i == 1){
          colnames(eem_df_ug_tofill) <- ex_scan_it_extrap$x
        }
        if(isTRUE(verbose)){
          message("finished loop ",i,"/",length(it_list))
        }
      }
    }
    eem_extrap <- eemdf_to_eem(eemdf = eem_df_ug_tofill,
                               file = eem$file,
                               sample = eem$sample,
                               location = eem$location) %>%
      'class<-'(c('eem'))
    return(eem_extrap)
  } else if(!is(eem,"eem")){
    stop("Please provide an object of class 'eem'")
  }
}
