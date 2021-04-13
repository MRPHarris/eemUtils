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
