#' Derive component Fmax values for a given PARAFAC model and eemlist.
#'
#' @description Fluorescence intensity at the component maxima (peak Ex/Em wavelength) can be used to infer intensity.
#'        Assuming a well-performing PARAFAC model, these values should correlate extremely well with the component
#'        A-mode values/loadings.
#'
#' @param pfmodel A single PARAFAC model object containing any number of components
#' @param eemlist a list of eems in the staRdom/eemR compliant format
#' @param component NULL or numeric. One or more components to extract fmax for. If NULL, all components targeted.
#' @param type two types of values are returned. 'fmax' for rescaled loadings (* BC mode maxima after Murphy et al., 2013), or "peakpick" for a per-sample intensity value picked at the BC mode maxima coordinates.
#' @param denormalise denormalise loadings prior to fmax calculation. Not necessary for peak-picking at the component spectra maxima.
#'
#' @export
#'
extrpf_fmax <- function(pfmodel, eemlist, component = NULL, type = "fmax", denormalise = FALSE){
  if(!is.null(component)){
    if(length(component) > 1){
      comps <- component
      # more than one component
      fmax_frame <- data.frame(matrix(NA,nrow = nrow(pfmodel$A), ncol = ncol(pfmodel$A[,c(component)])))
      colnames(fmax_frame) <- c(paste0("Comp.", comps))
    } else {
      # specific component chosen.
      comps <- component
      fmax_frame <- data.frame(matrix(NA,nrow = nrow(pfmodel$A), ncol = 1))
      colnames(fmax_frame) <- c(paste0("Comp.", component))
    }
  } else {
    # All components.
    comps <- seq(1,ncol(pfmodel$A),1)
    fmax_frame <- data.frame(matrix(NA,nrow = nrow(pfmodel$A), ncol = ncol(pfmodel$A)))
    colnames(fmax_frame) <- c(paste0("Comp.", comps))
  }
  fmax_frame$sample <- unlist(lapply(eemlist,"[[",'sample'))
  fmax_frame <- fmax_frame %>%
    dplyr::select("sample", everything())
  # By type.
  if(type == "fmax"){
    # denormalise, if need be
    if(isTRUE(denormalise)){
      loadings <- extrpf_loadings_denorm(pfmodel, eemlist)
    } else {
      loadings <- extrpf_loadings(pfmodel, eemlist)
    }
    # get the max B and C mode value
    maxvals <- vector("list", length = length(comps))
    for(cmp in seq_along(maxvals)){
      peak_it <- extrpf_peak_spectra(pfmodel, component = comps[cmp])
      mna <- (!is.na(peak_it$exn) & peak_it$ex == peak_it$max_ex[1])
      maxval_it <- peak_it[mna,]$value
      # now do this for each comp
      comp_it <- comps[cmp]
      compname_it <- paste0("Comp.",comp_it)
      loads_it <- loadings[,which(colnames(loadings) == compname_it)]
      fmax_frame[,which(colnames(fmax_frame) == paste0("Comp.",comp_it))] <- loads_it*maxval_it
    }
    return(fmax_frame)
  } else if(type == "peakpick"){
    # Get pfcomp peak positions
    complist <- vector("list", length = length(comps))
    peakpositions <- data.frame(matrix(NA,nrow = ncol(fmax_frame)), ncol = 2)
    colnames(peakpositions) <- c("Peak Excitation","Peak Emission")
    rownames(peakpositions) <- colnames(fmax_frame)
    for(c in seq_along(complist)){
      comp <- comps[c]
      spectra_it <- eemUtils::extrpf_peak_spectra(pfmodel, component = comp)
      peakpositions[c,2] <- spectra_it$max_em[1]
      peakpositions[c,1] <- spectra_it$max_ex[1]
    }
    # Extract maximum fluorescence value at peak maxima ex/em coordinates in each EEM
    for(f in seq_along(complist)){
      target_ex <- as.numeric(peakpositions[f,]$`Peak Excitation`)
      target_em <- as.numeric(peakpositions[f,]$`Peak Emission`)
      frame_it <- as.data.frame(matrix(NA,ncol = 1,nrow = length(eemlist)))
      colnames(frame_it) <- c("Intensity")
      for(e in seq_along(eemlist)){
        eem_it <- as.data.frame(eemlist[[e]], gather = FALSE)
        # index ex
        ex_vals <- as.numeric(colnames(eem_it))
        # index em
        em_vals <- as.numeric(rownames(eem_it))
        # get coords
        closest_ex <- as.numeric(which.min(abs(ex_vals - target_ex)))
        closest_em <- as.numeric(which.min(abs(em_vals - target_em)))
        # extract intensity value
        intensity_val_it <- eem_it[closest_em,closest_ex]
        # add it to frame
        frame_it[e,1] <- intensity_val_it
      }
      fmax_frame[,which(colnames(fmax_frame) == paste0("Comp.",comps[f]))] <- frame_it$Intensity
      message("Comp ",comps[f]," complete")
    }
    return(fmax_frame)
  } else {
    stop("Please supply type as either 'fmax' or 'peakpick'")
  }
}
