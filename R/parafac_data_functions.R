# Functions for data analysis or manipulation of PARAFAC model outputs.

#' Return a table of Core Consistency Diagnostic values for a given group of PARAFAC models
#'
#' @description A simple wrapper for stardom::eempf_corcondia().
#'
#' @param pfmodel A group of PARAFAC model objects, typically outputs from staRdom::eem_parafac()
#' @param eemlist The group of EEMs used to run the above PARAFAC model.
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
          "For information on the Core Consistency Diagnostic, see Bro and Kier (2003) and Murphy et al. (2013).")   # Completion message with refs
  if (any (Core_Consistencies <= 75)){
    message(" ------ \n","Warning: One or more of the models had a CCD <=75. This could indicate under-fitting of components.") # Warning message if CCD <=75
  }
  return(Core_Consistencies)
}

#' Extract fluorescence intensity loadings from a list of PARAFAC models.
#'
#' @description Extract the per-sample modeled intensity loadings for each component
#'       from a list of PARAFAC models. Optional numeric indexing of sample names
#'       using regmatching, assuming there are numerals to extract.
#'
#' @param pfmodel_list A list of PARAFAC models, typically an output from stardom::eem_parafac()
#' @param by_index TRUE/FALSE index the loadings by the sample names via registry matching. Requires numerals in the sample names.
#'
#' @export
#'
extrpf_loadings <- function(pfmodel_list, by_index = FALSE){
  loadings_list <- vector(mode = "list", length = length(pfmodel_list))
  pfmodel_name <- deparse(substitute(pfmodel_list))
  for (i in seq_along(pfmodel_list)){
    model <- pfmodel_list[[i]][["A"]]
    LoadingsA <- as.data.frame(model)
    Loadings <- as.data.frame(model)
    comp_num <- vector(mode = "list", length = ncol(Loadings))
    for (i in seq_along(comp_num)){
      cstat <- data.frame(matrix((paste0("C",i)),nrow(Loadings),1))
      Loadings <- cbind.data.frame(Loadings,cstat)
      col_index = i+(ncol(LoadingsA))
      colnames(Loadings)[col_index] <- paste0("C",i,"_Status")
    }
    data.table::setDT(Loadings, keep.rownames = TRUE)[]
    colnames(Loadings)[1] <- c("Sample_Name")
    if(isTRUE(by_index)){
      index <- as.numeric(regmatches(Loadings$Sample_Name, gregexpr("[[:digit:]]+", Loadings$Sample_Name)))
      Loadings <- cbind(Loadings, index)
    }
    loadings_list[[i]] <- Loadings
  }
  loadings_list <- loadings_list
  return(loadings_list)
}

#' Extract modeled spectra or EEM objects from a pfmodel list.
#'
#' @description Extracts either the component peak spectra, or the modelled EEMs,
#'     from a pfmodel list. Essentially produces the object used to plot the two
#'     types of staRdom::eempf_plot_comps().
#'
#' @param pfmodel_list PARAFAC model list.
#' @param type numeric 1 or 2 - extract EEMs or peaks, respectively.
#'
#' @export
#'
extrpf_spectra_or_eems <- function(pfmodel_list, type = 1){
  pfres <- pfmodel_list
  names = TRUE
  c <- pfres %>% lapply(eempf_comp_mat)
  colpal <- rainbow(75)[53:1]
  if (is.null(names(c))){
    names(c) <- paste0("model", seq(1:length(c)))
  }
  # Main results table, pre main rearrange
  tab <- lapply(1:length(c), function(n) {
    c1 <- c[[n]]
    mod_name <- names(c)[n]
    nc1 <- length(c1)
    nc2 <- 0
    lapply(c1, function(c2) {
      nc2 <<- nc2 + 1
      c2 <- c2 %>%
        mutate(comps = nc1,
               comp = paste0("Comp.", nc2),
               modname = mod_name)
    }) %>%
      bind_rows()
  }) %>%
    bind_rows() %>%
    mutate(modname = factor(modname, levels = names(c))) %>%
    mutate_at(vars(ex, em, value), as.numeric)
  fill_max <- tab$value %>% max(na.rm = TRUE)
  vals <- seq(from = 0, to = fill_max, length.out = length(colpal))
  vals <- (vals - min(vals))/diff(range(vals))
  if (type == 1) { # TYPE 1: Extract EEMs
    # Extract function would go here. Pull out items in a similar way to above.
    mod_num <- vector(mode = "list", length = length(pfmodel_list))
    mod_eem_list <- vector(mode = "list", length = length(pfmodel_list))
    model_names <- names(c)
    for (i in seq_along(mod_num)){
      Extracted_ModelEEMData <- tab %>%
        filter(comps == i & modname == (paste0("model",i)))
      mod_eem_list[[i]] = Extracted_ModelEEMData
      names(mod_eem_list)[i] = paste0(model_names[i],"_eems")
    }
    mod_eem_list <- mod_eem_list
  } else { # Type 2: Extract component maxima peak spectra
    plot_data <- tab %>%
      group_by(modname, comp) %>%
      mutate(max_pos = which.max(value), max_em = em[max_pos], max_ex = ex[max_pos]) %>%
      mutate(exn = ifelse(em == max_em, ex, NA), emn = ifelse(ex == max_ex, em, NA)) %>%
      filter(!is.na(emn) | !is.na(exn)) %>%
      ungroup()
    #plot_data <<- plot_data
    # Now to extract each set of data. Filter based upon number of components, I guess?
    mod_num <- vector(mode = "list", length = length(pfmodel_list))
    mod_spectra_list <- vector(mode = "list", length = length(pfmodel_list))
    model_names <- names(c)
    for (i in seq_along(mod_num)){
      Extracted_ModelSpectra <- plot_data %>%
        filter(comps == i & modname == (paste0("model",i)))
      mod_spectra_list[[i]] = Extracted_ModelSpectra
      names(mod_spectra_list)[i] = paste0(model_names[i],"_spectra")
    }
    mod_spectra_list <- mod_spectra_list
  }
}

#' Extract the modelled EEMs from a given PARAFAC model as an eemlist.
#'
#' @description Pulls out the modelled EEMs from a PARAFAC model, in the standard eemlist format,
#'      comprised of eem objects compliant with the staRdom/EEM/eemR framework.
#'
#' @param pfmodel a PARAFAC model object containing one or more components.
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
#' @description Determine the Ex/Em peak position of a PARAFAC component. Takes
#'      an output from extrpf_spectra_or_eems(type = 2).
#'
#' @param modeled_spectra PARAFAC-modeled spectra. An output from extrpf_spectra_or_eems(type = 2).
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

#' Derive the per-component percentage intensity contribution
#'
#' @description Determine the per-sample percentage contribution of each PARAFAC component
#'       across a given set of sample PARAFAC intensity loadings.
#'
#' @param loadings An output data frame from extrpf_loadings().
#'
#' @export
#'
get_pfload_percent <- function(loadings){
  ## n components
  ncomps <- as.numeric(length(unlist(which(grepl("Comp",colnames(loadings))))))
  ## which components are present, and in what columns do they lie?
  col_comps <- which(grepl("Comp",colnames(loadings)))
  what_comps <- unlist(lapply(str_split(colnames(loadings[col_comps]),"[.]"),"[[",2))
  # Colnames
  cnames_comps <- paste0("Percent_C",seq(1,ncomps,1))
  cnames <- c("Sample_Name",cnames_comps,"index")
  # Init empty frame
  pct_frame <- data.frame(matrix(nrow = nrow(loadings), ncol = length(cnames), data = NA))
  colnames(pct_frame) <- cnames

  # Assign sample names and indices
  pct_frame$Sample_Name <- loadings$Sample_Name
  pct_frame$index <- loadings$index

  # Total fluorescence
  fl_total <- rowSums(loadings[,col_comps])
  # Percent iteration loop
  it_list <- vector(mode = "list", length = ncomps)
  for(i in seq_along(it_list)){
    comp_it_pct <- (loadings[,col_comps[i]]/fl_total)*100
    pct_frame[,(1+i)] <- comp_it_pct
  }
  # Export data frame
  pct_frame
}

