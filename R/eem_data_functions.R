# These functions are related to the manipulation of EEM data.

#' Save a list of EEMs as individual .csv files.
#'
#' @description # Save the EEMs comprising an eemlist to a set of csv files. If no
#'      output_dir is set, the csvs are written straight to the working directory.
#'
#' @param eemlist A list of EEMs in a format compliant with eemR/staRdom.
#' @param outputfolder Either the path of the folder to write the .csvs to, or NULL for the working directory.
#' @param append_name optional; a character vector that will be appended to each .csv file.
#'
#' @export
#'
save_eemlist_csvs <- function(eemlist,outputfolder = NULL, append_name = NULL){
  if(is.null(outputfolder)){
    output_dir <- getwd()
    message("No outputfolder specified; exporting to working directory")
  }
  if(!is.null(append_name)){
    append_name = paste0("_",append_name)
  }
  for (i in seq_along(eemlist)){
    eem_ungathered <- as.data.frame(eemlist[[i]], gather = FALSE)
    write.csv(eem_ungathered, file = paste0(outputfolder,eemlist[[i]][["sample"]],append_name,".csv"))
    message(i,"/",length(eemlist)," | ",format(round((i/length(eemlist))*100,1), nsmall = 1),"% complete")
  }
  message("Done!")
}

#' Attempt to convert a raw EEM data dataframe to an eemR-compliant object of class 'eem'.
#'
#' @description Takes a data frame (e.g. raw EEM data) and attempts to coerce it into an eem object compliant with the eemR/staRdom framework.
#'
#' @param x The target dataframe.
#' @param sample Sample/EEM name.
#' @param location Optional; where is the EEM from
#' @param file Full file name and path of the EEM.
#' @param ex Which side is the excitation axis on - columns, or rows?
#'
#' @export
#'
data_frame_to_eem <- function(x, file = NULL, sample = NULL, location = NULL, ex = 'cols'){
  if(!(ex == 'cols' || ex == 'rows')){
    stop("ex must be either 'cols' or 'rows'")
  }
  eem <- vector(mode = 'list', length = 6)
  names(eem) <- c("file","sample","x","ex","em","location")
  # 6 parts to the list.
  eem[["file"]] <- file
  eem[["location"]] <- location
  eem[["sample"]] <- sample
  eem[["x"]] <- as.matrix(x)
  if(ex == 'cols'){
    eem[["ex"]] <- as.numeric(colnames(x))
    eem[["em"]] <- as.numeric(rownames(x))
  } else{
    eem[["ex"]] <- as.numeric(rownames(x))
    eem[["em"]] <- as.numeric(colnames(x))
  }
  e <- function(y){structure(y,class = 'eem')}
  eem <- e(eem)
}

#' Set all instances of '0'/zero to NA within a list of EEMs.
#'
#' @description Go through a list of eemR/staRdom-compliant EEMs, and replace all
#'       instances of 0 with NA. Used in cases where e.g. an Aqualog Rayleigh-masking
#'       step has masked those rows with 0, preventing proper interpolation and
#'       PARAFAC modeling.
#'
#' @param eemlist A list of EEMs in a format compliant with eemR/staRdom.
#' @param outputfolder Path to the folder where the new eemlist will be sent. Or NULL, if no export is desired.
#'
#' @export
#'
eem_0_to_NA <- function(eemlist, outputfolder = NULL){
  EEMs_DeNeg <- eemlist# load required packages
  for(i in seq_along(EEMs_DeNeg)){                                         # main for loop
    eem_ungathered <- as.data.frame(EEMs_DeNeg[[i]], gather = FALSE)       # extract EEM, don't gather
    eem_ungathered[,][eem_ungathered[,] == 0] <- NA                        # set all values of 0 in EEM to NA
    if(!is.null(outputfolder)){
      write.csv(eem_ungathered, file = paste0(outputfolder,EEMs_DeNeg[[i]][["sample"]],"_masked.csv"), row.names = TRUE) # Export EEM with iterative naming scheme.
    }
  }
  EEMs_DeNeg <- EEMs_DeNeg  # re-imports eems from the folder they were exported to
}

#' Set all negative values within a group of EEMs to 0.
#'
#' @description Negative fluorescence is not possible, and typically indicates
#'        bad processing, noise, or artefacts. This function will set all cases
#'        of negative fluorescence to 0.
#'
#' @param eemlist A list of EEMs in a format compliant with eemR/staRdom.
#' @param outputfolder Path to the folder where the new eemlist will be sent. Or NULL, if no export is desired.
#'
#' @export
#'
eem_neg_to_0 <- function(eemlist, outputfolder = NULL){
  EEMs_NoNeg <- eemlist
  for(i in seq_along(EEMs_NoNeg)){                                         # main for loop
    eem_ungathered <- as.data.frame(EEMs_NoNeg[[i]], gather = FALSE)       # extract EEM, don't gather
    eem_ungathered[,][eem_ungathered[,] <0] <- 0                        # set all values less than 0 in EEM to 0
    eemre <- data_frame_to_eem(eem_ungathered,
                      file = EEMs_NoNeg[[i]][['file']],
                      sample = EEMs_NoNeg[[i]][['sample']],
                      location = EEMs_NoNeg[[i]][['location']])
    EEMs_NoNeg[[i]] <- eemre
    if(!is.null(outputfolder)){
      write.csv(eem_ungathered, file = paste0(outputfolder,EEMs_NoNeg[[i]][["sample"]],"_noneg.csv"), row.names = TRUE) # Export EEM with iterative naming scheme.
    }
  }
  EEMs_NoNeg
}

#' Change wavelength range of eems with in an eemlist.
#'
#' @description A tweaked/re-built version of staRdom::eem_range(), which was not
#'      working for some datasets. Use this function to change the wavelength range
#'      of a group of EEMs.
#'
#' @param eemlist A list of EEMs compliant with the eemR/staRdom framework.
#' @param ex_range Numeric excitation min:max of the new range to constrict the EEMs to.
#' @param em_range Numeric excitation min:max of the new range to constrict the EEMs to.
#'
#' @export
#'
eem_range_mod <- function(eemlist, ex_range, em_range){
  if(class(eemlist) != "eemlist"){
    warning("please pass the function an object of class eemlist")
  }
  eemlist_name <- deparse(substitute(eemlist))
  # establish limits by indexing first EEM.
  # excitation
  min_ex <- as.numeric(min(ex_range)) #get min/max
  max_ex <- as.numeric(max(ex_range))
  ex_nm <- eemlist[[1]]$ex    # extract exitation wavelengths
  ex_nm <- ex_nm[length(ex_nm):1]   # flip to allow indexing
  dt_ex <- data.table(ex_nm, val = ex_nm)   # data.table
  setattr(dt_ex,"sorted", "ex_nm")   # set data.table attributes so it can be indexed
  min_ex_index_r <- as.numeric(dt_ex[data.table(max_ex), roll = "nearest", which = TRUE]) # get min index
  min_ex_index <- length(eemlist[[1]]$ex) - (min_ex_index_r - 1)   # flip, as the true wavelengths are the other way round
  max_ex_index_r <- as.numeric(dt_ex[data.table(min_ex), roll = "nearest", which = TRUE]) # same as above, for max
  max_ex_index <- length(eemlist[[1]]$ex) - (max_ex_index_r - 1)
  # emission
  min_em <- as.numeric(min(em_range))
  max_em <- as.numeric(max(em_range))
  em_nm <- eemlist[[1]]$em   # acquiring nearest value in table.
  dt_em <- data.table(em_nm, val = em_nm)
  setattr(dt_em, "sorted", "em_nm")
  min_em_index <- as.numeric(dt_em[data.table(min_em), roll = "nearest", which = TRUE]) # indem of said value
  max_em_index <- as.numeric(dt_em[data.table(max_em), roll = "nearest", which = TRUE])
  is_between <- function(x, a, b) {
    x >= a & x <= b
  }
  new_eemlist <- eemlist
  for (i in seq_along(eemlist)){
    # excitation range
    if (!missing(ex_range)){
      stopifnot(is.numeric(ex_range), all(ex_range >= 0))
      index_ex <- which(is_between(eemlist[[i]]$ex, eemlist[[i]]$ex[max_ex_index],
                                   eemlist[[i]]$ex[min_ex_index]))
      if (length(index_ex != 0)) {
        eemlist[[i]]$ex <- eemlist[[i]]$ex[index_ex]
        eemlist[[i]]$x <- eemlist[[i]]$x[, index_ex]
      }
    }
    # emission range
    if (!missing(em_range)){
      stopifnot(is.numeric(em_range), all(em_range >= 0))
      index_em <- which(is_between(eemlist[[i]]$em, eemlist[[i]]$em[min_em_index],
                                   eemlist[[i]]$em[max_em_index]))
      if (length(index_em != 0)) {
        eemlist[[i]]$em <- eemlist[[i]]$em[index_em]
        eemlist[[i]]$x <- eemlist[[i]]$x[index_em, ]
      }
    }
    new_eemlist[[i]] <- eemlist[[i]]
  }
  new_eemlist
}

#' Set areas of one or more EEMs to NA.
#'
#' @description A tweaked/re-built version of staRdom::set_NA(). Sets areas of
#'        eems to NA. Option to interpolate. Changes from set_NA() relate to how
#'        the range is identified. I recommend interpolating seperately.
#'
#' @param eemlist A list of EEMs in a format compliant with eemR/staRdom.
#' @param sample optional, names or indices of samples to process
#' @param em optional, emission wavelengths to set NA
#' @param ex optional, excitation wavelengths to set NA
#' @param interpolate FALSE, 1 or 2, interpolate NAs or not, 2 different methods.
#'
#' @import data.table
#'
#' @export
#'
eem_setNA_mod <- function(eem_list, sample = NULL, em = NULL, ex = NULL, interpolate = FALSE){
  if (is.null(sample)) {
    sample <- eem_names(eem_list)
  }
  if (is.numeric(sample)) {
    sample <- eem_names(eem_list)[sample]
  }
  eem_list <- lapply(eem_list, function(eem) {
    if (eem$sample %in% sample) {
      if (is.null(ex)){
        ex2 <- 1:ncol(eem$x)
      } else (ex2 <- which(eem$ex %in% ex))
      if (is.null(em)){
        em2 <- 1:nrow(eem$x)
      }
      if(length(em) > 1){
        min_em <- as.numeric(min(em))
        max_em <- as.numeric(max(em))
        em_nm <- eem$em   # acquiring nearest value in table.
        dt_em <- data.table(em_nm, val = em_nm)
        setattr(dt_em, "sorted", "em_nm")
        min_em_val <- as.numeric(dt_em[J(min_em), roll = "nearest"][1,2]) # Identify value closest to that specified
        min_em_index <- as.numeric(dt_em[J(min_em), roll = "nearest", which = TRUE]) # index of said value
        max_em_val <- as.numeric(dt_em[J(max_em), roll = "nearest"][1,2])
        max_em_index <- as.numeric(dt_em[J(max_em), roll = "nearest", which = TRUE])
        em <- eem$em[min_em_index:max_em_index]
        em2 <- which(eem$em %in% em)}  # new em comprised of actual values
      if(length(em) == 1){
        em_nm <- eem$em
        dt_em <- data.table(em_nm, val = em_nm)
        setattr(dt_em, "sorted", "em_nm")
        em_val <- as.numeric(dt_em[J(em), roll = "nearest"][1,2])
        em_index <- as.numeric(dt_em[J(em), roll = "nearest", which = TRUE])
        em2 <- which(eem$em %in% em)}
      eem$x[em2, ex2] <- NA
    }
    eem
  }) %>% `class<-`("eemlist")
  if (interpolate != FALSE) {
    eem_list[which(eem_names(eem_list) %in% sample)] <- eem_interp(eem_list[which(eem_names(eem_list) %in%
                                                                                    sample)], type = interpolate) %>% `class<-`("eemlist")
  }
  eem_list
}

#' Normalise EEM data using one of a number of different methods.
#'
#' @description Normalise EEMs in a dataset using different methods: convert from
#'      raman peak to raman area (Raman Units), using values supplied from a list,
#'      using a single value, or internally normalise each EEM by its own maximum intensity.
#'
#' @param eemlist A list of EEMs, compliant with the eemR/staRdom framework.
#' @param type one of 'raman_peak_to_area','standard_from_list','standard_single_value', or 'max_eem_val'.
#' @param norm_log A dataframe listing all eems in the eemlist, along with which sample group they belong to.
#' @param norm_data A dataframe listing the normalisation values (old and new) corresponding to each sample group.
#' @param norm_value If a single value is being used, a numeric value. If not, leave as NULL.
#' @param output_folder Path for export of normalised EEMs.
#' @param append_name Extension to be added to the name of each EEM.
#'
#' @export
#'
normalise_eemlist <- function(eemlist, type = 'raman_peak_to_area', norm_log, norm_data, norm_value = NULL,
                              outputfolder = NULL, append_name = NULL){
  # Available types:
  # - raman_peak_to_area: denormalise by area of the RamanPeak, renormalise by the peak area. Most must be supplied.
  # - standard_from_list: supplied a single value, via a norm_values dataframe.
  # - standard_single_value: a single value that all the EEMs will be normalised by.
  # - max_eem_val: normlises each EEM in the eemlist by the highest intensity of each individual EEM.

  # This function requires 3 input data sets (if using any type other than standard_single_value):
  # 1. an eemlist
  # 2. dataframe (norm_log) listing all eems (samples) in the eemlist, along with which sample group they belong to.
  # 3. dataframe (norm_values) listing the normalisation values (old and new) corresponding to each sample group.
  # input checks
  if(!(type == 'raman_peak_to_area' || type == 'standard_from_list' || type == 'standard_single_value' || type == 'max_eem_val')){
    stop("Please specify type as one of 'raman_peak_to_area','standard_from_list','standard_single_value','max_eem_val'")
  }
  if(!(class(eemlist) == 'eemlist')){ # make sure the eemlist is, in fact, an eemlist
    stop("The input eemlist must be a list of class 'eemlist'")
  }
  if(!is.null(append_name)){ # Check name appending, and set up future eemlist append
    if(!(class(append_name) == "character")){
      stop("If supplying a name to append to output files, please ensure it is a valid character string (i.e. of class 'character')")
    }
    append_name_forlist = paste0("_",append_name)
  } else {
    append_name_forlist = "_normalised"
  }
  #required functions
  # removed data_frame_to_eem() here, as it's now a part of the package.
  # removed save_eemlist_csvs(), as it's part of the package.
  # setup
  normalised_eemlist <- vector(mode = 'list',length = length(eemlist))
  e <- function(y){structure(y,class = 'eemlist')}
  normalised_eemlist <- e(normalised_eemlist)
  # type loops. Deal with exporting separately.
  if(type == 'raman_peak_to_area'){
    for(i in seq_along(eemlist)){
      # extract sample name, filename, location
      sample_name <- eemlist[[i]][["sample"]]
      filename <- eemlist[[i]][["file"]]
      location <- eemlist[[i]][["location"]]
      # ungather
      eem_ungathered <- as.data.frame(eemlist[[i]], gather = FALSE)
      # Extract relevant rows from norm_log and norm_data frame
      target_norm_log <- norm_log[grep(sample_name, norm_log$sample_name),]
      sample_group <- as.numeric(target_norm_log$sample_group) # identify which sample group the sample belongs to
      target_norm_data_row <- as.numeric(which(norm_data$sample_group %in% sample_group))
      target_norm_data <- norm_data[target_norm_data_row,]
      # get the norm values
      raman_peak_area <- target_norm_data[,2]
      raman_peak_emission <- target_norm_data[,3]
      # do the denormalising and normalising
      eem_denorm <- eem_ungathered*raman_peak_emission
      eem_norm <- eem_denorm/raman_peak_area
      # coerce back to an eem object
      eem_norm <- data_frame_to_eem(x = eem_norm, sample = sample_name, file = filename, location = location, ex = 'cols')
      # insert into normalised eemlist
      normalised_eemlist[[i]] <- eem_norm
      #write.csv(eem_ungathered, file = paste0(outputfolder,eemlist[[i]][["sample"]],"_denorm.csv"), row.names = TRUE)
      message(i,"/",length(eemlist)," complete | type = raman_peak_to_area, norm_value = ", raman_peak_area)
    }
  }
  if(type == 'standard_from_list'){
    for(i in seq_along(eemlist)){
      sample_name <- eemlist[[i]][["sample"]]
      filename <- eemlist[[i]][["file"]]
      location <- eemlist[[i]][["location"]]
      # ungather
      eem_ungathered <- as.data.frame(eemlist[[i]], gather = FALSE)
      # Extract relevant rows from norm_log and norm_data frame
      target_norm_log <- norm_log[grep(sample_name, norm_log$sample_name),]
      sample_group <- as.numeric(target_norm_log$sample_group) # identify which sample group the sample belongs to
      target_norm_data_row <- as.numeric(which(norm_data$sample_group %in% sample_group))
      target_norm_data <- norm_data[target_norm_data_row,]
      # get the norm values
      norm_value <- target_norm_data[,2]
      # do the denormalising and normalising
      eem_norm <- eem_ungathered/norm_value
      # coerce back to an eem object
      eem_norm <- data_frame_to_eem(x = eem_norm, sample = sample_name, filename = filename, location = location, ex = 'cols')
      # insert into normalised eemlist
      normalised_eemlist[[i]] <- eem_norm
      #write.csv(eem_ungathered, file = paste0(outputfolder,eemlist[[i]][["sample"]],"_denorm.csv"), row.names = TRUE)
      message(i,"/",length(eemlist)," complete | type = standard_from_list, norm_value = ", norm_value)
    }
  }
  if(type == 'standard_single_value'){
    if(is.null(norm_value)){ # check that a norm value has been specified
      stop("If using a single value normalisation, please supply the norm_value as an input parameter.")
    }
    for(i in seq_along(eemlist)){ # perform normalisation on each eem within the eemlist
      # extract sample name, filename, location
      sample_name <- eemlist[[i]][["sample"]]
      filename <- eemlist[[i]][["file"]]
      location <- eemlist[[i]][["location"]]
      # ungather, conduct normalisation
      eem_ungathered <- as.data.frame(eemlist[[i]], gather = FALSE)
      eem_norm <- eem_ungathered/norm_value
      # coerce back to an eem object
      eem_norm <- data_frame_to_eem(x = eem_norm, sample = sample_name, filename = filename, location = location, ex = 'cols')
      # insert into normalised eemlist
      normalised_eemlist[[i]] <- eem_norm
      #write.csv(eem_ungathered, file = paste0(outputfolder,eemlist[[i]][["sample"]],"_denorm.csv"), row.names = TRUE)
      message(i,"/",length(eemlist)," complete | type = standard_single_value, norm_value = ",norm_value,"")
    }
  }
  if(type == 'max_eem_val'){
    for(i in seq_along(eemlist)){
      sample_name <- eemlist[[i]][["sample"]]
      filename <- eemlist[[i]][["file"]]
      location <- eemlist[[i]][["location"]]
      # ungather
      eem_ungathered <- as.data.frame(eemlist[[i]], gather = FALSE)
      # Extract relevant rows from norm_log and norm_data frame
      # get the norm values
      norm_value <- max(eem_ungathered, na.rm = TRUE)
      # do the denormalising and normalising
      eem_norm <- eem_ungathered/norm_value
      # coerce back to an eem object
      eem_norm <- data_frame_to_eem(x = eem_norm, sample = sample_name, file = filename, location = location, ex = 'cols')
      # insert into normalised eemlist
      normalised_eemlist[[i]] <- eem_norm
      #write.csv(eem_ungathered, file = paste0(outputfolder,eemlist[[i]][["sample"]],"_denorm.csv"), row.names = TRUE)
      message(i,"/",length(eemlist)," complete | type = max_eem_val, norm_value = ", norm_value)
    }
  }
  # Output checks
  normalised_eemlist <- plyr::compact(normalised_eemlist)
  # Exporting, if a directory is specified.
  if(!is.null(outputfolder)){
    message("Exporting EEMs to output folder")
    save_eemlist_csvs(normalised_eemlist, outputfolder = outputfolder, append_name = append_name) #save eemlist
  }
  # Outputs
  return(normalised_eemlist)
  message("Normalisation complete")
}

#' Average a set of EEMs.
#'
#' @description Average the EEMs within an eemlist.
#'
#' @param eemlist A list of EEMs, compliant with the eemR/staRdom framework.
#'
#' @export
#'
eemlist_average <- function(eemlist){
  if(length(eemlist) == 1){
    message("1 EEM passed to average_eems() for averaging. Returning unchanged.")
    new_eemlist <- eemlist
    return(new_eemlist)
  } else if(length(eemlist) > 1){
    ungathered_list <- vector(mode = "list",length = length(eemlist))
    for(i in seq_along(eemlist)){
      eem_it <- eemlist[[i]]
      file_it <- eem_it[['file']]
      sample_it <- eem_it[['sample']]
      location_it <- eem_it[['location']]
      eem_ungathered <- as.data.frame(eemlist[[i]], gather = FALSE)       # extract EEM, don't gather
      ungathered_list[[i]] <- eem_ungathered
    }
    # Average it
    averaged <- Reduce("+",ungathered_list)/length(ungathered_list)
    # Now back to eem
    averaged_eem <- eemdf_to_eem(averaged,
                                 file = "",
                                 sample = "averaged_eem",
                                 location = "")
    new_eemlist <- vector(mode = "list",length = 1)
    class(new_eemlist) <- "eemlist"
    message(paste0(length(eemlist)," EEMs passed to average_eems() for averaging."))
    new_eemlist[[1]] <- averaged_eem
    return(new_eemlist)
  }
}

#' Takes a data frame and attempts to coerce it to an EEM object of the style used by EEM/eemR/staRdom.
#'
#' @description An alternative to eemR's eem constructor. Mirrored across from the SampleQueue package utility functions. Required for eemlist_average().
#'
#' @param eemdf the dataframe to be coerced to an EEM object.
#' @param file filename of the EEM, if applicable.
#' @param sample the samplename of the EEM, if applicable.
#' @param location the location of the EEM file, if applicable.
#'
#' @noRd
eemdf_to_eem <- function(eemdf,
                         file,
                         sample,
                         location){
  # code adapted from staRdom's .eem_csv importer.
  x <- eemdf
  ex <- colnames(x)[] %>% as.numeric()
  em <- rownames(x) %>% as.numeric()
  x <- x[,] %>% as.matrix() %>% unname()
  x <- x[!is.na(em),!is.na(ex)]
  ex <- ex[!is.na(ex)]
  em <- em[!is.na(em)]
  l <- list(
    file = file,
    sample = sample,
    x = x,
    ex = ex,
    em = em,
    location = location
  )
  class(l) <- "eem"
  return(l)
}
