# Functions for data analysis or manipulation of PARAFAC model outputs.

#' Return a table of Core Consistency Diagnostic values for a given group of PARAFAC models
#'
#' @description A simple wrapper for stardom::eempf_corcondia().
#'
#' @param pfmodel A group of PARAFAC model objects, typically outputs from staRdom::eem_parafac()
#' @param eemlist The group of EEMs used to run the above PARAFAC model
#'
#' @export
#'
Generate_CORCONDIA <- function(pfmodel,eemlist){
  message(length(pfmodel)," models detected. Spinning up...") # Initialization message
  Core_Consistencies <- data.frame(matrix(0,length(pfmodel)))             # blank frame for CCD values
  row_names_temp <- data.frame(matrix(0,length(pfmodel)))                 # blank frame for row names, to be indexed to # components/model
  for(i in seq_along(pfmodel)){                                           # initiate main loop
    CC_values <- staRdom::eempf_corcondia(pfmodel[[i]],eemlist)                      # Calculate iterative CCD for each model in list
    Core_Consistencies <- rbind.data.frame(Core_Consistencies, CC_values)      # bind CCD stats to frame
    row_names_iter <- ncol(pfmodel[[i]][["A"]])                                      # get # cols in model loading frame (proxy for # components)
    row_names_temp <- rbind.data.frame(row_names_temp, row_names_iter)         # bind row name (# components) to frame
    Sys.sleep(0.01)                                                            # pause ever so briefly
    if (CC_values <= 75){                                                      # warning threshold set to CC of <= 75. Tweak to suit!
      message("Heads up: CCD for the ",row_names_iter,"-component model is <=75!") # Warning message if CCD <=75
    }
    cat(" ------ \n",i,"/",length(pfmodel), "Complete","\n","------ \n")             # print iteration complete/total iterations
  }
  Core_Consistencies <- Core_Consistencies[-(1:length(pfmodel)), , drop = FALSE]                   # Drop blank rows
  row_names_temp <- row_names_temp[-(1:length(pfmodel)), , drop = FALSE]                           # Drop blank rows
  row.names(Core_Consistencies) <- paste0(row_names_temp[[1]]," Component/s")                # apply indexed row names to the CCD frame
  colnames(Core_Consistencies) <- c('Core Consistency')                                      # column name for the CCD frame
  print(Core_Consistencies)                                                                  # print the CCD frame
  message("\n","Core Consistency Diagnostic/s generated for all models within ", deparse(substitute(pfmodel)),".","\n",
          "For information on the Core Consistency Diagnostic (and its limitations), see Bro and Kier (2003) and Murphy et al. (2013).")   # Completion message with refs
  if (any (Core_Consistencies <= 75)){
    message(" ------ \n","Warning: One or more of the models had a CCD <=75. This could indicate over-fitting of components.") # Warning message if CCD <=75
  }
  return(Core_Consistencies)
}

#' Extract fluorescence intensity loadings from a PARAFAC model.
#'
#' @description Extract the per-sample modeled intensity loadings for each component
#'       from a PARAFAC models. Optional numeric indexing of sample names
#'       using regmatching, assuming there are numerals to extract. Use lapply for
#'       multiple models.
#'
#' @param pfmodel A PARAFAC model object containing any number of components
#' @param by_index TRUE/FALSE index the loadings by the sample names via registry matching. Requires numerals in the sample names
#'
#' @export
#'
extrpf_loadings <- function(pfmodel, by_index = FALSE, fill_stat = FALSE){
  model <- pfmodel[["A"]]
  LoadingsA <- as.data.frame(model)
  Loadings <- as.data.frame(model)
  comp_num <- vector(mode = "list", length = ncol(Loadings))
  if(isTRUE(fill_stat)){
    for (i in seq_along(comp_num)) {
      cstat <- data.frame(matrix((paste0("C", i)),nrow(Loadings), 1))
      Loadings <- cbind.data.frame(Loadings, cstat)
      col_index = i + (ncol(LoadingsA))
      colnames(Loadings)[col_index] <- paste0("C",i, "_Status")
    }
  }
  # data.table::setDT(Loadings, keep.rownames = TRUE)[]
  Loadings <- Loadings %>%
    rownames_to_column("sample") %>%
    select(sample, everything())
  colnames(Loadings)[1] <- c("sample")
  if (isTRUE(by_index)) {
    index <- as.numeric(regmatches(Loadings$sample,gregexpr("[[:digit:]]+", Loadings$sample)))
    Loadings <- cbind(Loadings, index)
  }
  Loadings
}

#' Multiply the loadings values of PARAFAC components by normalisation factors
#'
#' @description A PARAFAC model fed with normalised eems outputs by default loadings that are less useful for direct
#'       sample-to-sample comparison. This can be remedied by multiplying the loadings by each eem fmax value (i.e. the
#'       normalisation factors used to normalise the eems originally).
#'
#' @param pfmodel a PARAFAC model object, typically an individual ouptut from staRdom::eem_parafaC()
#' @param eemlist a list of eem objects compliant with the staRdom/eemR framework
#' @param type short or long format. Short by default. Long format data works better for grouping in ggplot
#'
#' @export
#'
extrpf_loadings_denorm <- function(pfmodel, eemlist, type = "short"){
  maxvals <- eemlist_fmax_values(eemlist)
  fm <- extrpf_loadings(pfmodel)
  loadings_frame <- fm[,2:ncol(fm)]
  newframe <- apply(loadings_frame, 2, function(col) {
    col * maxvals
  }) %>% data.frame()
  newframe <- newframe %>% `colnames<-`(c(paste0("Comp.",
                                                 seq(1, ncol(newframe), 1)))) %>% `rownames<-`(unlist(lapply(eemlist,
                                                                                                             "[[", "sample"))) %>% rownames_to_column(var = "sample")
  if (type == "short") {
    newframe
  }
  else if (type == "long") {
    newframe <- newframe %>% pivot_longer(cols = c(2:ncol(newframe)))
  }
  else {
    stop("Unknown 'type'. Please input either 'short' or 'long'")
  }
}


#' Extract the modelled EEMs from a given PARAFAC model as an eemlist.
#'
#' @description Pulls out the modelled EEMs from a PARAFAC model, in the standard eemlist format,
#'      comprised of eem objects compliant with the staRdom/EEM/eemR framework.
#'
#' @param pfmodel a PARAFAC model object containing one or more components
#'
#' @export
#'
extrpf_eems <- function(pfmodel){
  # Code based upon staRdom::eempft_comp_mat()
  eemlist <- lapply(seq(1:ncol(pfmodel$A)), function(comp){
    m <- matrix(pfmodel$B[, comp]) %*% t(matrix(pfmodel$C[,comp])) %>% data.frame()
    colnames(m) <- pfmodel$C %>% rownames()
    rownames(m) <- pfmodel$B %>% rownames()
    m_eem <- eemUtils::eemdf_to_eem(m, file = NULL, sample = colnames(pfmodel$B)[comp], location = NULL, gathered = FALSE)
    m_eem
  })
  names(eemlist) <- colnames(pfmodel$B)
  eemlist
}

#' Determines the peak positions for a PARAFAC modeled component.
#'
#' @description DEFUNCT Determine the Ex/Em peak position of a PARAFAC component. Takes
#'      an output from extrpf_spectra_or_eems(type = 2).
#'
#' @param modeled_spectra PARAFAC-modeled spectra. An output from extrpf_spectra_or_eems(type = 2)
#'
#' @export
#'
peak_positions <- function(modeled_spectra){
  peak_list <- vector(mode = "list", length = length(modeled_spectra))
  for(i in seq_along(modeled_spectra)){
    spectra <- modeled_spectra[[i]]
    ncomponents <- as.numeric(unique(spectra[["comps"]]))
    peak_coords <- data.frame(matrix(NA,ncomponents,2))
    maxem <- rle(spectra[["max_em"]])[[2]]
    maxex <- rle(spectra[["max_ex"]])[[2]]
    peak_coords[1:length(maxex),1] = maxex
    peak_coords[1:length(maxem),2] = maxem
    rownames_vals <- seq(1:(ncomponents))
    rownames(peak_coords) = paste("Component_",rownames_vals,sep = "")
    col1 = paste0("Peak Excitation")
    col2 = paste0("Peak Emission")
    colnames(peak_coords) = c(col1, col2)
    # add to list instead
    peak_list[[i]] <- peak_coords
  }
  print(peak_list)
  peak_list <- peak_list
}


#' Extract spectra at the peak wavelength position.
#'
#' @description An alternative to the extrpf_peaks_or_eems. Simply extract the spectra at the peak for a given component.
#'
#' @param pfmodel A PARAFAC model object a la staRdom
#' @param component Extract spectra for this component. Numeric
#'
#' @export
#'
extrpf_peak_spectra <- function(pfmodel, component = 1){
  pfres <- pfmodel # Code adapted from staRdom package; naming will reflect this.
  names = TRUE
  c <- eempf_comp_mat(pfres) # get matrices of all the comps.
  if (is.null(names(c))) {
    names(c) <- paste0("model", seq(1:length(c)))
  }
  # Performing for only one model.
  c1 <- c # holdover code; lapply wrap removed here.
  nc1 <- length(c1)
  nc2 <- 0
  # for each component, pull out the data
  tab <- lapply(c1, function(c2) {
    nc2 <<- nc2 + 1
    c2 <- c2 %>% mutate(comps = nc1, comp = paste0(nc2))
  }) %>% bind_rows() %>%
    mutate_at(vars(ex, em, value, comp), as.numeric) # numeric vars
  # collate, organise, extract max spectra.
  comp_spectra <- tab %>%
    group_by(comp) %>%
    mutate(max_pos = which.max(value), max_em = em[max_pos], max_ex = ex[max_pos]) %>%
    mutate(exn = ifelse(em == max_em, ex, NA), emn = ifelse(ex == max_ex, em, NA)) %>%
    dplyr::filter(!is.na(emn) | !is.na(exn)) %>%
    ungroup() %>%
    mutate_at(vars(ex, em, value, comp), as.numeric) %>%
    dplyr::filter(comp == component)
  comp_spectra
}

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
    select("sample", everything())
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

#' Extract the per-EEM percentage contributions of each component
#'
#' @description Return the percent contribution of each component to the modelled fluorescence response of each sample, using
#'      the modelled loadings
#'
#' @param pfmodel a PARAFAC model object
#' @param eemlist supply the eemlist used for denormalisation. Only necessary if denormalise is set to TRUE.
#' @param denormalise TRUE/FALSE to denormalise the loadings prior to the percent calculation
#'
#' @export
#'
extrpf_loadings_percent <- function(pfmodel, eemlist, denormalise = FALSE){
  if(isTRUE(denormalise)){
    loadings <- extrpf_loadings_denorm(pfmodel, eemlist)
  } else {
    loadings <- extrpf_loadings(pfmodel)
  }
  FI_totals <- rowSums(loadings[,2:ncol(loadings)])
  pct_contrib <- loadings[,2:ncol(loadings)] %>%
    mutate_all(.,function(col){(col/FI_totals)*100})
  pct_contrib$sample <- loadings$sample
  pct_contrib <- pct_contrib %>%
    select('sample',everything())
  pct_contrib
}

