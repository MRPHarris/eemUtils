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

# Normalise all data within an eemlist. Options as 'types'. Expand function for info.
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
      eem_norm <- data_frame_to_eem(x = eem_norm, sample = sample_name, filename = filename, location = location, ex = 'cols')
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
      eem_norm <- data_frame_to_eem(x = eem_norm, sample = sample_name, filename = filename, location = location, ex = 'cols')
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
  normalised_eemlist <<- normalised_eemlist
  assign(paste0(deparse(substitute(eemlist)), append_name_forlist), normalised_eemlist,  envir = parent.frame())  # new name for eemlist based upon input eemlist
  rm(normalised_eemlist, envir = parent.frame())  # remove repeat data
  message("Normalisation complete")
}

