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
}

