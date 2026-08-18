#' Extract spectra from an EEM at a given Ex/Em wavelength.
#'
#' @description 'Slice' an EEM at a given Ex/Em coordinate pair, and extract the spectra at that location.
#'
#' @param eem an eem object compliant with the staRdom/eemR framework
#' @param ex an excitation wavelength value, in nm
#' @param em an emission wavelength value, in nm
#'
#' @export
#'
slice_eem <- function(eem, ex, em){
  # input checks
  if (is.data.frame(eem)) {
  } else if (class(eem) == "eem") {
    eem_df <- as.data.frame(eem, gather = FALSE)
  } else {
    stop("Please pass the function an object of class 'eem'")
  }
  # Extract emission values slice at given ex wavelength
  if(ex %in% colnames(eem_df)){
    em_slice <- data.frame(matrix(NA, nrow = nrow(eem_df), ncol = 2))
    em_slice[, 1] <- as.numeric(rownames(eem_df))
    em_slice[, 2] <- as.numeric(as.matrix(eem_df[, which(colnames(eem_df) == ex)]))
    colnames(em_slice) <- c("emission", "intensity")
    em_slice <- pivot_longer(em_slice, cols = emission, values_to = "wavelength")
  } else {
    ex <- binary_search_nearest(data = colnames(eem_df), value = ex)
    em_slice <- data.frame(matrix(NA, nrow = nrow(eem_df), ncol = 2))
    em_slice[, 1] <- as.numeric(rownames(eem_df))
    em_slice[, 2] <- as.numeric(as.matrix(eem_df[, which(colnames(eem_df) == ex)]))
    colnames(em_slice) <- c("emission", "intensity")
    em_slice <- pivot_longer(em_slice, cols = emission, values_to = "wavelength")
  }
  # Extract excitation values slice at given em wavelength
  if(em %in% rownames(eem_df)){
    # Continue as normally
    ex_slice <- data.frame(matrix(NA, nrow = ncol(eem_df), ncol = 2))
    ex_slice[, 1] <- as.numeric(colnames(eem_df))
    ex_slice[, 2] <- as.numeric(as.matrix(t(eem_df[which(rownames(eem_df) == em), ])))
    colnames(ex_slice) <- c("excitation", "intensity")
    ex_slice <- pivot_longer(ex_slice, cols = excitation, values_to = "wavelength")
  } else {
    # Oh no! Find nearest em value with binary search via data.table.
    em <- binary_search_nearest(data = rownames(eem_df), value = em)
    ex_slice <- data.frame(matrix(NA, nrow = ncol(eem_df), ncol = 2))
    ex_slice[, 1] <- as.numeric(colnames(eem_df))
    ex_slice[, 2] <- as.numeric(as.matrix(t(eem_df[which(rownames(eem_df) == em), ])))
    colnames(ex_slice) <- c("excitation", "intensity")
    ex_slice <- pivot_longer(ex_slice, cols = excitation, values_to = "wavelength")
  }
  # Combine
  slices <- rbind(em_slice, ex_slice)
  slices %>% mutate_at(vars(intensity, wavelength), as.numeric)
  slices
}
