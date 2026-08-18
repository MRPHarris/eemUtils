#' Import and parse a downloaded export .csv file from the OpenFluor databse
#'
#' @description The OpenFluor database allows for the statistical comparison of a PARAFAC model with published spectra.
#'   These comparisons can be downloaded as a .csv file, which this function can parse into a sorted list object containing all
#'   the information within the file. It is recommended to first open and immediately save the downloaded .csv in e.g. MS Excel as the
#'   raw file is slightly finnicky with read.csv and read_csv.
#'
#' @param file full filepath to the .csv file.
#' @param normalise_spectra TRUE/FALSE; logical to normalise excitation and emission spectra from the database (irrelevant/automatic during the TCC comparisons but can help with plotting)
#'
#' @importFrom readr read_csv
#' @importFrom dplyr relocate
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#'
#' @export
#'
openfluor_parse <- function(file, normalise_spectra = TRUE){
  ## read file
  file <- read_csv(file, col_names = FALSE, guess_max = 1000)
  ## Is the datset name too long?
  if(!is.na(nchar(file[2,2]) > 0)){
    warning('The OpenFluor dataset name is longer than 70 characters, which may cause problems when parsing spectra. This can be avoided by shortening the dataset name in OpenFluor prior to comparisons and export.')
  }
  ## List object
  parse_list <- vector('list', length = 8)
  names(parse_list) <- c('Dataset / Model Name','Export Date',
                         'TCCex',"TCCem","TCCexem",
                         'Matches','Match info', 'Match spectra')
  ### Sections
  ## Information
  nchars <- function(x){
    x2 <- mapply(nchar,x)
    x2
  }
  initial_lines <- split(file[c(2:3),1], seq(nrow(file[c(2:3),1])))
  parse_list[c(1,2)] <- lapply(initial_lines, function(line){
    # line <- initial_lines[[2]]
    line2 <- line[which(line != "NA")] %>% unlist()
    line3 <- paste(line2[which(nchars(line2) > 0)], collapse = "")
    line4 <- trimws(paste(unlist(lapply(strsplit(line3 ,"[:]"),"[",c(2:length(unlist(strsplit(line3 ,"[:]")))))), collapse = ":"), which = 'both')
    line4
  })
  ## Match criteria
  first_sumline <- as.character(file[as.numeric(which(mapply(grepl,"Summary",file[,1]))[1]),1])
  rmv <- c("[(]","[)]")
  parse_list[[3]] <- as.numeric(str_remove_all(trimws(unlist(lapply(strsplit(unlist(lapply(strsplit(first_sumline,"TCCex"),"[[",2)),"[-]"),"[[",1))),paste(rmv,collapse = "|")))
  parse_list[[4]] <- as.numeric(str_remove_all(trimws(unlist(lapply(strsplit(unlist(lapply(strsplit(first_sumline,"TCCem"),"[[",2)),"[-]"),"[[",1))),paste(rmv,collapse = "|")))
  parse_list[[5]] <- parse_list[[3]] * parse_list[[4]]
  ## Matching models
  # Identify which lines read 'models"
  model_list_start <- as.numeric(which(file[,1] == 'Models'))
  model_list_end <- as.numeric(which(mapply(grepl,"Summary",file[,1]))[1])-2
  models <- data.frame(matrix(NA,nrow = (model_list_end - model_list_start), ncol = 2)) %>% 'colnames<-'(c("Entry","Identifier"))
  models$Entry <- file[c(seq(model_list_start+1,model_list_end,1)),1]
  models$Identifier <- file[c(seq(model_list_start+1,model_list_end,1)),2]
  parse_list[[6]] <- models
  ## Model match information
  # Which rows start with "Summary of"?
  summary_headers <- as.numeric(which(mapply(grepl,"Summary",file[,1])))
  if(length(summary_headers) != nrow(models)){
    warning("Mismatch between summary entries and matching models. Some data may not have been imported from OpenFluor; check .csv file.")
  }
  matchedcomp_headers <- as.numeric(which(mapply(grepl,"Matched Component",file[,1])))
  nsumdatrows <- matchedcomp_headers[1] - summary_headers[1]
  nmatchedcomps <- (matchedcomp_headers[1] - summary_headers[1]) - length(summary_headers)*3
  match_info <- data.frame(matrix(NA,nrow = nmatchedcomps, ncol = 6)) %>%
    'colnames<-'(c('matching_dataset','comp_internal','comp_external','TCCex','TCCem','TCCexem'))
  # Iterate over entries, adding them to match_info
  match_entries <- vector('list',length = length(summary_headers))
  for(m in seq_along(match_entries)){
    # m = 1
    # message(paste0('starting iteration ',m))
    header <- as.character(file[summary_headers[m],1])
    info <- as.character(file[summary_headers[m]+1,1])
    if(m != length(match_entries)){
      match_it <- file[c(seq(summary_headers[m]+2,summary_headers[m+1]-2,1)),]
    } else {
      match_it <- file[c(seq(summary_headers[m]+2,matchedcomp_headers[1]-3,1)),]

    }
    # Trim NA
    match_trimmed <- match_it[,colSums(is.na(match_it))<nrow(match_it)]
    # Trim whitespace
    no_ws_check <- apply(match_trimmed,2,function(col){
      # col <- match_trimmed[,6]
      colcheck <- which(col == '')
      if(length(colcheck) != 0){
        r <- FALSE
      } else {
        r <- TRUE
      }
      r
    })
    match_ws_trimmed <- match_trimmed[,no_ws_check]
    match_ws_trimmed$entry <- trimws(unlist(lapply(strsplit(unlist(lapply(strsplit(header,"Matches for "),"[[",2)),"[-]"),"[[",1)))
    match_ws_trimmed <- match_ws_trimmed %>%
      relocate(entry)
    match_entries[[m]] <- match_ws_trimmed
  }
  match_entries <- rlist::list.rbind(match_entries) %>%
    'colnames<-'(c('matching_dataset','comp_internal','comp_external','TCCex','TCCem','TCCexem'))
  parse_list[[7]] <- match_entries

  ## Matched component spectra
  matched_comp_spectra <- vector('list', length = length(matchedcomp_headers))
  names(matched_comp_spectra) <- trimws(unlist(lapply(strsplit(c(file[matchedcomp_headers,1])[[1]],"[:]"),"[[",2)))
  # iterate over each component, formatting excitation and emission spectra
  for(comp in seq_along(matched_comp_spectra)){
    comp = 9
    matched_headerline <- as.character(file[matchedcomp_headers[comp],1])
    if(comp != length(matched_comp_spectra)){
      complines_it <- matchedcomp_headers[comp + 1] - matchedcomp_headers[comp]
      complines <- file[c(seq(matchedcomp_headers[comp]+2,matchedcomp_headers[comp]+complines_it-2,1)),]
    } else {
      complines_it <- nrow(file) - matchedcomp_headers[comp]
      complines <- file[c(seq(matchedcomp_headers[comp]+2,matchedcomp_headers[comp]+complines_it,1)),]
    }
    ex_n <- as.numeric(which(mapply(grepl,"Type: Ex",complines[,1])))
    em_n <- as.numeric(which(mapply(grepl,"Type: Em",complines[,1])))
    # Excitation spectra
    ex_lines <- complines[c(seq(ex_n+1,em_n-1,1)),]
    ex_lines <- replace(ex_lines,ex_lines=="",NA)
    em_lines <- complines[c(seq(em_n+1,nrow(complines),1)),]
    em_lines <- replace(em_lines,em_lines=="",NA)
    # Now parse all spectra.
    ## Excitation
    ex_names <- ex_lines[,1]
    ex_dat <- ex_lines %>%
      select(-c(X1)) %>%
      mutate(across(where(is.factor), as.character)) %>%
      mutate(across(where(is.character),as.numeric)) %>%
      t() %>%
      'colnames<-'(c(ex_names)[[1]])%>%
      data.frame()
    em_names <- em_lines[,1]
    em_dat <- em_lines %>%
      select(-c(X1)) %>%
      mutate(across(where(is.factor), as.character)) %>%
      mutate(across(where(is.character),as.numeric)) %>%
      t() %>%
      'colnames<-'(c(em_names)[[1]]) %>%
      data.frame()
    if(isTRUE(normalise_spectra)){
      em_wl <- em_dat[,1]
      em_dat <- apply(em_dat[2:ncol(em_dat)],2, function(col){
        col <- col/max(col, na.rm = T)
      }) %>% data.frame() %>% mutate(wavelength = em_wl) %>% relocate(wavelength)
      ex_wl <- ex_dat[,1]
      ex_dat <- apply(ex_dat[2:ncol(ex_dat)],2, function(col){
        col <- col/max(col, na.rm = T)
      }) %>% data.frame() %>% mutate(wavelength = ex_wl) %>% relocate(wavelength)
    }
    complist <- vector('list', length = 2) %>% 'names<-'(c('excitation','emission'))
    complist[[1]] <- em_dat
    complist[[2]] <- ex_dat
    matched_comp_spectra[[comp]] <- complist
  }
  parse_list[[8]] <- matched_comp_spectra
  return(parse_list)
}
