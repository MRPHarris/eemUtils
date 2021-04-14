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
  IR = zeros(NoSamples,2)
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
      p = zeros(length(sg),1)
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
