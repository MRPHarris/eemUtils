# Various miscellaneous functions for eem or fluorometer datasets.
# Here are all the functions that didn't fit in elsewhere.

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

#' Extract the peak areas from spectrometer Raman curves, for EEM normalisation purposes.
#'
#' @description This function provides two methods to calculate the area under the
#'       Raman peak of water for fluorescence data normalisation purposes. The first
#'       method is a port of RamanIntegrationRange from the MATLAB drEEM package,
#'       which uses a gradient detection method to identify the start and end of
#'       the Raman peak of water. The second is a standard Aqualog method, which
#'       simply assumes that the peak extends from 380 to 410nm. Plot export and
#'       .gif wrapping are optional. Some example raman curves are included in
#'       the eemUtils data folder.
#'
#' @param RAMdat A data frame containing wavelength increments in column 1, and fluorescence intensities in subsequent columns.
#' @param range_lower optional; a lower limit for calculations.
#' @param range_upper optional; an upper limit for calculations. 500 by default to avoid higher-order scattering bands.
#' @param method one of 'aqualog' or 'RIR' for either default Aqualog 380:410nm, or to use RamanIntegrationRange's gradient calculations to determine peak position, and thus area.
#' @param halfwidth if method = 'RIR', what is the halfwidth? 1800 by default.
#' @param sequence if method = 'RIR', what is the sequence? 8 by default.
#' @param tolerance if method = 'RIR', what is the gradient tolerance? 0.05 by default.
#' @param landa what is the excitation band of the raman spectra? 350 (nm) by default.
#' @param output_dir optional; path to the output directory of peak images.
#' @param gif optional; create .gif animation of peak images?
#'
#' @export
#'
extract_ramanpeak_areas <- function(RAMdat, range_lower = NULL, range_upper = 500, method = "aqualog",halfwidth = 1800, sequence = 8,
                                    tolerance = 0.05, landa = 350, output_dir = NULL, gif = FALSE){
  ### SETUP
  landa_ex <- landa
  # Input checks
  if(!(method == "aqualog" || method == "RIR")){
    warning("method must be one of 'aqualog' or 'RIR'. Defaulting to 380:410 IR aqualog method.")
    method = "aqualog"
  }
  if(range_upper > 650){
    message("Warning: an emission range extending past ~650nm may cause the function to select for the secondary Rayleigh","\n",
            "scatter peak, which occurs at ~700nm emission at 350nm excitation. Set the range_upper parameter","\n",
            "to <= 650 to avoid this.")
  }
  if(is.null(range_lower)){
    range_lower <- RAMdat[1,1]
  }
  if(is.null(range_upper)){
    range_upper <- RAMdat[1,ncol(RAMmat)]
  }
  embound_lower <- which.min(abs(RAMdat[,1] - range_lower)) # lower bounds
  embound_upper <- which.min(abs(RAMdat[,1] - range_upper)) # upper bounds
  RAMdat <- RAMdat[embound_lower:embound_upper,] # new WL range applied.
  new_min_ex <- round(min(RAMdat[,1], na.rm = TRUE))
  new_max_ex <- round(max(RAMdat[,1], na.rm = TRUE))
  new_WL_range <- unique(seq(new_min_ex,new_max_ex,1)) # new range with distance of 1
  check_diffs <- diff(new_WL_range) # acquire running diffs
  ### SLIT WIDTH INTERPOLATION
  # Interpolating across wavelength range to ensure slit width is 1 going forward
  new_em_check <- approx(x = RAMdat[,1], # linear interpolation using approx()
                         y = RAMdat[,2],
                         xout = new_WL_range,
                         method = "linear")
  new_em_check_df <- data.frame(matrix(t(unlist(new_em_check)), nrow = length(new_em_check$x), ncol = 2))
  interp_df <- data.frame(matrix(NA,nrow = nrow(new_em_check_df), ncol = ncol(RAMdat)))
  interp_df[,1] <- new_WL_range
  colnames(interp_df)[1] <- "wavelength"
  nSamp <- ncol(RAMdat)-1
  qit_list <- vector(mode = "list", length = nSamp)
  for(i in seq_along(qit_list)){
    new_em <- approx(x = RAMdat[,1], y = RAMdat[,i+1], xout = new_WL_range, method = "linear") # interpolate emission values across new wavelength range for sample i
    new_em <- data.frame(matrix(t(unlist(new_em)), nrow = length(new_em$x), ncol = 2)) # unlist to dataframe
    interp_df[,i+1] <- new_em[,2] # add interp data to interp_df
    colnames(interp_df)[i+1] <- c(paste0("sample_",i)) # column name to label sample (where sample i is in column i+1)
  } # interp loop
  # From this point onward, the interpolated values will be used instead, and will be transposed to rows. NAs removed.
  RAMmat <- interp_df
  RAMmat <- t(RAMmat)
  col_ok = apply(RAMmat,2,function(x)!any(is.na(x))) # identify columns containing NAs, assign logical T/F
  RAMmat <- RAMmat[,col_ok] # remove NAs
  #plot(x = RAMmat[1,], y = RAMmat[2,])
  # confirm transposing. If not, transpose!
  if(length(RAMmat) < nrow(RAMmat)){ # check that the data is in rows. If not, transpose.
    RAMmat = t(RAMmat)
  } # Ensure R is in rows, not columns. Transposes if in columns.
  # Use ispos to confirm wavelengths are in the first row of RAMmat.
  ispos <- function(Xmatrix){
    has.neg <- apply(as.data.frame(Xmatrix), 1, function(col) any(col <= 0))
    if(any(has.neg == TRUE, na.rm = TRUE)){
      TrueOrFalse = 0
    } else {
      TrueOrFalse = 1
    }
  }
  isneg <- function(Xmatrix){
    has.pos <- apply(as.data.frame(Xmatrix), 1, function(col) any(col > 0))
    if(any(has.pos == TRUE, na.rm = TRUE)){
      TrueOrFalse = 0
    } else {
      TrueOrFalse = 1
    }
  }
  if(ispos(Xmatrix = pracma::gradient(RAMmat[1,]))==0){
    warning("check that wavelengths are in the first row of RAMmat, and that scans for different samples are in subsequent rows!")
  }
  R <- RAMmat
  # Some themes for later plots
  ggplot2::theme_set(cowplot::theme_cowplot(12))
  geom.text.size = 3
  theme.size = (14/5) * geom.text.size
  inset.geom.text.size = 2
  inset.theme.size = (14/5) * inset.geom.text.size
  theme_paper_text <- theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = theme.size),
    axis.text = element_text(size = theme.size),
    axis.text.x = element_text(),
    axis.text.y = element_text(
      hjust = 0.8,
      margin = margin(r = 5)),
    axis.ticks.length = unit(.25, "cm"),
    legend.background = element_blank(),
    legend.key.height = unit(7,"pt"),
    legend.spacing.y = unit(7, "pt"),
    legend.text = element_text(size = theme.size),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA))
  # Number of samples
  NoSamples = as.numeric(length(RAMmat[,1])) - 1 # n samples
  it_list <- vector(mode = "list", length = NoSamples) # iteration list for main loop
  IR = matrix(0, nrow = NoSamples, ncol = 2)
  # The following code was originally wrapped in the 'RamanPeakPos' function from drEEM. Please cite accordingly.
  landa_em = 1/(1/landa_ex-3380/10^7)
  landa_em <<- landa_em
  Emmin = 1/(1/landa_em + (halfwidth*10^-7))
  Emmax = 1/(1/landa_em - (halfwidth*10^-7))
  Emmin = floor(Emmin)
  Emmax = ceiling(Emmax)
  Emmin = Emmin
  Emmax = Emmax
  message("Legacy peak max at em",round(landa_em*10)/10,", lying between ex",Emmin,":",Emmax)
  landa_em = round(landa_em*10)/10
  # end RamanPeakPos wrap
  RU_Norm_Fac_df <- data.frame(matrix(NA,nrow = NoSamples, ncol = 2))
  RU_Norm_Fac_df[,1] <- seq(1,NoSamples,1)
  colnames(RU_Norm_Fac_df) <- c("sample_group","raman_peak_area")
  # commence main loop over rows.
  if(method == "aqualog"){
    for(z in seq_along(it_list)){
      zn <- formatC(z, width = 2, format = "d", flag = "0") # formatted sample numbers (for <100 sample groups)
      R <- RAMmat[c(1,z+1),] # pull out sample data.
      x = R[1,]
      y = R[2,]
      IR[z,1] <- 380
      IR[z,2] <- 410
      bounds_A <- which.min(abs(R[1,] - IR[z,1]))
      bounds_B <- which.min(abs(R[1,] - IR[z,2]))
      old_bounds_A <- as.numeric(which.min(abs(R[1,] - Emmin)))
      old_bounds_B<- as.numeric(which.min(abs(R[1,] - Emmax)))
      # Integration Range
      integration_spectra <- as.data.frame(t(R))
      integration_spectra <- as.data.frame(integration_spectra[bounds_A:bounds_B,])
      colnames(integration_spectra)[2] <- "values"
      #plot(R[1,old_bounds_A:old_bounds_B],R[2,old_bounds_A:old_bounds_B], type = 'l')
      legacy_spectra <- as.data.frame(t(R))
      legacy_spectra <- as.data.frame(legacy_spectra[old_bounds_A:old_bounds_B,])
      colnames(legacy_spectra)[2] <- "values"
      # Baseline trapezoid. from start of legacy peak.
      #trapezoid <- legacy_spectra[c(1,nrow(legacy_spectra)),]
      #trap_slope <- diff(trapezoid$values)/diff(trapezoid$wavelength)
      #x <- trapezoid[,1]
      #y <- trapezoid[,2]
      #trap <- approx(x = x, y = y, xout = seq(min(x),max(x),1))
      #trap_df <- data.frame(matrix(unlist(trap), nrow=length(trap$x), ncol = length(trap)))
      #colnames(trap_df) <- c("x","y")
      # Baseline trapezoid. From start of aqualog peak.
      trapezoid <- integration_spectra[c(1,nrow(integration_spectra)),]
      trap_slope <- diff(trapezoid$values)/diff(trapezoid$wavelength)
      x <- trapezoid[,1]
      y <- trapezoid[,2]
      trap <- approx(x = x, y = y, xout = seq(min(x),max(x),1))
      trap_df <- data.frame(matrix(unlist(trap), nrow=length(trap$x), ncol = length(trap)))
      colnames(trap_df) <- c("x","y")
      # trapezoid bounds that intersect with the integrated area.
      trap_bounds_A <- which.min(abs(trap_df[,1] - IR[z,1]))
      trap_bounds_B <- which.min(abs(trap_df[,1] - IR[z,2]))
      BL_df <- trap_df[trap_bounds_A:trap_bounds_B,]
      colnames(BL_df) <- c("x","y")
      # Raman Peak Area
      RU_norm_factor <- sum(integration_spectra[,2]) # sum intensities under the peak.
      RU_norm_factor <- RU_norm_factor - sum(BL_df$y) # subtract baseline.
      RUnorm_disp <- round(RU_norm_factor,5)
      RU_Norm_Fac_df[z,2] <- RU_norm_factor
      message("[",zn,"] IR = ",IR[z,1],":",IR[z,2], " | Peak Area = ",RU_norm_factor) #integration range. PRINTS
      # plot
      if(!(is.null(output_dir))){
        png(paste0(output_dir,"Sampleset",zn,"_aqualog.png"), units="cm", width=10, height=10, res=300)
        plot <- ggplot() +
          geom_line(data = legacy_spectra, aes(x = wavelength, y = values)) +
          geom_bar(data = integration_spectra, aes(x = wavelength, y = values), stat = 'identity', fill = "#828282",colour = "#262626", width = 0.7, size = 0.2) +
          geom_line(data = trap_df, aes(x = x, y = y), colour = 'red', size = 1) +
          labs(title = paste0("Sample set ",z),
               x = "Wavelength (nm)",
               y = "Calibrated intensity (Sc/Rc)") +
          scale_y_continuous(limits = c(-5,200), expand = c(0,0)) +
          theme_paper_text +
          theme(plot.title = element_text(hjust = 0.5)) +
          annotate("text", x = 405, y = 150, hjust = 0, size = inset.geom.text.size, label = paste0("raman peak area = ",RUnorm_disp)) +
          annotate("text", x = 405, y = 142, hjust = 0, size = inset.geom.text.size, label = paste0("tolerance = ",tolerance*100,"%")) +
          annotate("text", x = 405, y = 136, hjust = 0, size = inset.geom.text.size, label = deparse(bquote(halfwidth == .(halfwidth)~cm^-1)), parse = T) +
          annotate("text", x = IR[z,1]-1, y = 30, hjust = 0, size = inset.geom.text.size, label = IR[z,1]) +
          annotate("text", x = IR[z,2], y = 30, hjust = 0, size = inset.geom.text.size, label = IR[z,2])
        print(plot)
        dev.off()
        print(plot)
      }
    }
  }
  if(method == "RIR"){
    for(z in seq_along(it_list)){
      zn <- formatC(z, width = 2, format = "d", flag = "0") # formatted sample numbers (for <100 sample groups)
      tolcheck <- tolerance
      R <- RAMmat[c(1,z+1),] # pull out sample data.
      x = R[1,]
      y = R[2,]
      g <- gradient(R[2,]) # determine gradient, return y output.
      sg = forecast::ma(g, 5, TRUE) #5 point moving average
      sg2 = smooth(g)
      sg <- as.matrix(sg)
      #Identify sequences of positive and negative gradients
      p = matrix(0,nrow=length(sg),ncol=1)
      n = p
      sq = sequence - 1
      for (i in seq_along(1:length(sg))){
        sgsub_pos = sg[i]
        p[i] = ispos(sgsub_pos)
      } # logical 1/0 for positives
      for (i in seq_along(1:length(sg))){
        sgsub_neg = sg[i]
        n[i] = isneg(sgsub_neg)
      } # logical 1/0 for negatives
      sgp = sg*p
      sgn = sg*n
      # restrict gradient bounds to avoid rayleigh scatter in determining the maximum gradient
      sgmax = max(sg[x > Emmin], na.rm = TRUE)
      # incorporate tolerance
      xi = x[sgp - tolcheck*sgmax > 0]
      sgpi = sgp[sgp - tolcheck*sgmax > 0]
      xf = x[sgn + tolcheck*sgmax < 0]
      sgnf = sgn[sgn + tolcheck*sgmax < 0]
      #positive sequences within integration range (before peak)
      gradpos <- matrix(NA, nrow = 2, ncol = length(xi[xi > Emmin & xi < landa_em]))
      gradpos[1,] <- xi[xi > Emmin & xi < landa_em]
      gradpos[2,] <- sgpi[xi > Emmin & xi < landa_em]
      gradpos <- as.data.frame(t(gradpos))
      pos_diffs <- diff(gradpos[,1])# acquire diffs
      pos_diffs[length(pos_diffs)+1] <- pos_diffs[length(pos_diffs)]
      if(length(which(pos_diffs > 1))){ # remove those rows from gradpos
        pos_diff_index <- as.numeric(which(pos_diffs > 1)+1)
        gradpos <- gradpos[max(pos_diff_index):nrow(gradpos),]
      }
      gradpos <- as.matrix(t(gradpos))
      col_ok_pos = apply(gradpos,2,function(x)!any(is.na(x))) # identify columns containing NAs, assign logical T/F
      gradpos <- gradpos[,col_ok_pos]
      #plot(gradpos[1,],gradpos[2,])
      # Negative sequences within integration range (after peak).
      gradneg <- matrix(NA, nrow = 2, ncol = length(xf[xf > landa_em & xf < Emmax]))
      gradneg[1,] <- xf[xf > landa_em & xf < Emmax]
      gradneg[2,] <- sgnf[xf > landa_em & xf < Emmax]
      gradneg <- as.data.frame(t(gradneg))
      neg_diffs <- diff(gradneg[,1])# acquire diffs for values > peak position. Gaps in the positive gradient are excluded.
      neg_diffs[length(neg_diffs)+1] <- neg_diffs[length(neg_diffs)]
      if(length(which(neg_diffs > 1))){ # remove those rows from gradpos
        neg_diff_index <- as.numeric(which(neg_diffs > 1))
        gradneg <- gradneg[1:min(neg_diff_index),]
      }
      gradneg <- as.matrix(t(gradneg))
      col_ok_neg = apply(gradneg,2,function(x)!any(is.na(x))) # identify columns containing NAs, assign logical T/F
      gradneg <- gradneg[,col_ok_neg]
      #plot(gradneg[1,],gradneg[2,])
      # Messages for lack of pos/neg values.
      if(is.null(gradpos)){
        message("Heads up: No positive sequences of desired length within the integration range.","\n",
                "          Check that excitation wavelengths are correctly specified.","\n",
                "          To trouble shoot this error, try reducing sequence length and setting tolerance to 0.")
      }
      if(is.null(gradneg)){
        message("Heads up: No negative sequences of desired length within the integration range.","\n",
                "          Check that excitation wavelengths are correctly specified.","\n",
                "          To trouble shoot this error, try reducing sequence length and setting tolerance to 0.")
      }
      # Compose Integration Range
      IR[z,1] <- floor(gradpos[1,1]*2)/2
      IR[z,2] <- ceiling(gradneg[1,ncol(gradneg)]*2)/2
      bounds_A <- which.min(abs(R[1,] - IR[z,1]))
      bounds_B <- which.min(abs(R[1,] - IR[z,2]))
      #plot(R[1,bounds_A:bounds_B],R[2,bounds_A:bounds_B], type = 'l') # plot the defined Raman curve
      # The 'old peak' as defined by RamanPeakPos().
      old_bounds_A <- as.numeric(which.min(abs(R[1,] - Emmin)))
      old_bounds_B<- as.numeric(which.min(abs(R[1,] - Emmax)))
      #plot(R[1,old_bounds_A:old_bounds_B],R[2,old_bounds_A:old_bounds_B], type = 'l')
      legacy_spectra <- as.data.frame(t(R))
      legacy_spectra <- as.data.frame(legacy_spectra[old_bounds_A:old_bounds_B,])
      gradient_spectra <- as.data.frame(t(R))
      gradient_spectra <- as.data.frame(gradient_spectra[bounds_A:bounds_B,])
      colnames(legacy_spectra)[2] <- "values"
      colnames(gradient_spectra)[2] <- "values"
      #plot(legacy_spectra$wavelength,legacy_spectra$values, type = 'l') # plot 370:428 (legacy peak)
      #plot(gradient_spectra$wavelength,gradient_spectra$values, type = 'l')
      # calculate baseline trapezoid
      trapezoid <- gradient_spectra[c(1,nrow(gradient_spectra)),]
      trap_slope <- diff(trapezoid$values)/diff(trapezoid$wavelength)
      x <- trapezoid[,1]
      y <- trapezoid[,2]
      BL <- approx(x = x, y = y, xout = seq(min(x),max(x),1))
      BL_df <- data.frame(matrix(unlist(BL), nrow=length(BL$x), ncol = length(BL)))
      colnames(BL_df) <- c("x","y")
      # RU normalisation factor
      RU_norm_factor <- sum(gradient_spectra[,2]) # sum intensities under the peak.
      RU_norm_factor <- RU_norm_factor - sum(BL_df$y) # subtract baseline.
      RUnorm_disp <- round(RU_norm_factor,5)
      RU_Norm_Fac_df[z,2] <- RU_norm_factor
      message("[",zn,"] IR = ",IR[z,1],":",IR[z,2], " | Peak Area = ",RU_norm_factor) #integration range. PRINTS
      # plot
      if(!(is.null(output_dir))){
        png(paste0(output_dir,"Sampleset",zn,"_RIR.png"), units="cm", width=10, height=10, res=300)
        plot <- ggplot() +
          geom_line(data = legacy_spectra, aes(x = wavelength, y = values)) +
          geom_bar(data = gradient_spectra, aes(x = wavelength, y = values), stat = 'identity', fill = "#828282",colour = "#262626", width = 0.7, size = 0.2) +
          geom_line(data = BL_df, aes(x = x, y = y), colour = 'red', size = 1) +
          labs(title = paste0("Sample set ",z),
               x = "Wavelength (nm)",
               y = "Calibrated intensity (Sc/Rc)") +
          scale_y_continuous(limits = c(-5,200), expand = c(0,0)) +
          theme_paper_text +
          theme(plot.title = element_text(hjust = 0.5)) +
          annotate("text", x = 405, y = 150, hjust = 0, size = inset.geom.text.size, label = paste0("raman peak area = ",RUnorm_disp)) +
          annotate("text", x = 405, y = 142, hjust = 0, size = inset.geom.text.size, label = paste0("tolerance = ",tolerance*100,"%")) +
          annotate("text", x = 405, y = 136, hjust = 0, size = inset.geom.text.size, label = deparse(bquote(halfwidth == .(halfwidth)~cm^-1)), parse = T) +
          annotate("text", x = IR[z,1]-1, y = 30, hjust = 0, size = inset.geom.text.size, label = IR[z,1]) +
          annotate("text", x = IR[z,2], y = 30, hjust = 0, size = inset.geom.text.size, label = IR[z,2])
        print(plot)
        dev.off()
        print(plot)
      }
    }
  }
  integration_ranges <<- IR # output IR
  raman_peak_areas <<- RU_Norm_Fac_df # output RU normalisation factors
  assign(paste0("raman_peak_areas_",method), raman_peak_areas,  envir = parent.frame())
  rm(raman_peak_areas, envir = parent.frame())
  if(isTRUE(gif)){ # create gif
    if(!is.null(output_dir)){
      list.files(path = output_dir, pattern = paste0("*",method,".png"), full.names = TRUE) %>%
        image_read() %>% # reads each path file
        image_join() %>% # joins image
        image_animate(fps=4) %>% # animates, can opt for number of loops
        image_write(paste0(output_dir,"Peaks_",method,".gif"))
    } else{
      message("No output_dir specified; gif not exported.")
    }
  }
}

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

#' Get SSC along with alpha and beta penalty terms
#'
#' @description The shift- and shape sensitive congruence (SSC) was developed by Wunsch et al., 2019 as an improvement
#'      upon the TCC metric. It incorporates two penalty terms, alpha and beta, to account for differences in the wavelength
#'      peak position and area. This function adds these terms to the data frame returned by staRdom::ssc()
#'
#' @param mat1 a matrix
#' @param mat2 a matrix
#' @param tcc TRUE/FALSE to return only TCC value instead of SSC, alpha and beta.
#'
#' @export
#'
ssc_more <- function (mat1, mat2, tcc = FALSE) {
  if (any(is.null(mat1), is.na(mat1), is.null(mat2), is.na(mat2))) {
    a <- NA
  } else {
    a <- lapply(1:ncol(mat1), function(nc) {
      col1 <- mat1[, nc]
      apply(mat2, 2, function(col2) {
        tcc_cal <- sum(col1 * col2)/sqrt(sum(col1^2) *
                                           sum(col2^2))
        if (!tcc) {
          wl <- as.numeric(names(col1))
          if (any(is.na(wl)) | pracma::isempty(wl)) {
            stop("SSCs cannot be calculated. Please add wavelengths as rownames of the matrices!")
          }
          alpha <- abs((wl[which.max(col1)] - wl[which.max(col2)])/diff(range(wl)))
          beta <- abs((sum(col1/max(col1)) - sum(col2/max(col2)))/diff(range(wl)))
          ssc <- tcc_cal - alpha - beta
          ssc <- c(ssc, alpha, beta)
          # rownames(a) <- (c("ssc","alpha","beta"))
        } else {
          tcc_cal
        }
      })
    }) %>% setNames(colnames(mat1)) %>% do.call(rbind, .)
  }
  attr(a, "method") <- ifelse(tcc, "TCC", "SSC")
  if(!isTRUE(tcc)){
    rownames (a) <- c("ssc","alpha","beta")
  }
  a
}

#' Plot two spectra and show congruence metric data
#'
#' @description plots two matrices (e.g. spectra) in a similar fashion to figure 2 from Wunsch et al., 2019.
#'
#' @param mat1 a matrix
#' @param mat2 a matrix
#' @param normalise TRUE/FALSE to normalise the matrices before comparison and subsequent plotting.
#'
ssc_plot <- function(mat1, mat2, normalise = TRUE){
  mat_labels <- c(deparse(substitute(mat1)), deparse(substitute(mat2)))
  if(!all(names(mat1) == names(mat2))){
    message("Warning: there is a misalignment in the matrix wavelengths.")
  }
  if(isTRUE(normalise)){
    mat1 <- mat1/max(mat1)
    mat2 <- mat2/max(mat2)
  }
  # Get SSC and associated vars
  ssc_fp <- ssc_more(mat1 = mat1, mat2 = mat2, tcc = FALSE)
  tcc_fp <- ssc(mat1 = mat1, mat2 = mat2, tcc = TRUE)
  # format data for plotting.
  SSC_table <- data.frame(matrix(cbind(mat1,mat2,rownames(mat1), rownames(mat2)), ncol = 4))
  colnames(SSC_table) <- c("mat1_value","mat2_value","mat1_wavelength", "mat2_wavelength")
  SSC_table <- SSC_table %>%
    pivot_longer(everything(), names_sep = "_", names_to = c("mat",".value")) %>%
    mutate_at(vars(value,wavelength), as.numeric)
  # colours
  clrs_fill <- c("#6b6b6b","#3e80a8")
  clrs_line <- c("#383838","#234961")
  # plotting
  plt <- ggplot() +
    geom_line(data = SSC_table, aes(x = wavelength, y = value, group = mat, colour = mat, linetype = mat), size = 0.8) +
    theme_cowplot(12) +
    geom_area(data = SSC_table, aes(x = wavelength, y = value, fill = mat, group = mat), alpha = 0.5, position = "identity") +
    scale_linetype_manual(values = c(1, 2), labels = mat_labels) +
    scale_color_manual(values = clrs_line, labels = mat_labels) +
    scale_fill_manual(values = clrs_fill, labels = mat_labels) +
    scale_y_continuous(expand = c(0,0), limits = c(0,(max(SSC_table$value)+(max(SSC_table$value)/10)))) +
    scale_x_continuous(expand = c(0,0)) +
    labs(fill = "Spectra", color = "Spectra", linetype = "Spectra",
         x = "Wavelength") +
    annotate(geom = "text", x = (max(SSC_table$wavelength)-(max(SSC_table$wavelength)/10)), y = (max(SSC_table$value)-(max(SSC_table$value)/10)),
             label = paste0("TCC = ",    format(round(tcc_fp,3), nsmall = 3), "\n",
                            "SSC = ",    format(round(ssc_fp[1,],3), nsmall = 3), "\n",
                            "\u03B1 = ", format(round(ssc_fp[2,],3), nsmall = 3),"\n",
                            "\u03B2 = ", format(round(ssc_fp[3,],3), nsmall = 3)),
             hjust = 1)
  if(isTRUE(normalise)){
    plt <- plt +
      labs(y = "Normalised value")
  }
  plt
}

