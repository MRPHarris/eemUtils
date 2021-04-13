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
