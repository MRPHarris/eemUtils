#' An importer for Aqualog-produced 2D emission scans
#'
#' @description A simple importing function for 2D emission scans produced using the
#'      HORIBA Aqualog. Supports both sample and blank scan types (though not both at the same time)
#'
#' @param files character vector containing one or more full filenames of emission scans
#' @param type character - one of either 'sample' or 'blank'
#'
#' @export
#'
emscan_read <- function (files, type = "sample"){
  if (!(type == "sample") && !(type == "blank")) {
    stop("Please supply param 'type' as either 'blank' or 'sample'.")
  }
  if (!is.list(files)) {
    files <- as.list(files)
  }
  if (type == "sample") {
    emscans <- lapply(files, function(f) {
      fname <- unlist(lapply(str_split(f, "[/]"),
                             tail, n = 1L))
      data <- readLines(f)
      dat_read <- read.delim(f)
      emscan <- vector("list", length = ncol(dat_read)+2)
      names(emscan) <- c("name","type",colnames(dat_read))
      names(emscan)[grep('Sample..',names(emscan))] <- 'emission'
      for(i in seq_along(dat_read)){
        col_it <- dat_read[3:nrow(dat_read),i]
        vals <- as.numeric(col_it[which(!is.na(col_it))])
        vals <- as.numeric(vals[which(!is.na(vals))])
        emscan[[i+2]] <- vals
      }
      emscan[["name"]] <- fname
      emscan[["type"]] <- "sample"
      emscan
    })
  } # Adding a note
  else if (type == "blank") {
    emscans <- lapply(files, function(f) {
      fname <- unlist(lapply(str_split(f, "[/]"),
                             tail, n = 1L))
      data <- readLines(f)
      dat <- stringr::str_extract_all(data, "-?\\d+(?:\\.\\d*)?(?:[eE][+\\-]?\\d+)?")
      headline <- unlist(str_split(data[1], "[\t]"))
      headline2 <- unlist(str_split(data[3], "[\t]"))
      headline_combined <- c(headline[1], headline2[2],
                             headline2[3], headline[4], headline[5], "scan_ex",
                             headline2[7], headline2[8], headline[9], headline[10],
                             headline2[11], headline[12], headline2[13], "corrected emission")
      emscan <- as.list(as.numeric(unlist(dat[4])))
      names(emscan) <- headline_combined
      emscan[[1]] <- as.numeric(unlist(lapply(dat[4:length(dat)],
                                              "[[", 1)))
      emscan[[2]] <- as.numeric(unlist(lapply(dat[4:length(dat)],
                                              "[[", 2)))
      emscan[[3]] <- as.numeric(unlist(lapply(dat[4:length(dat)],
                                              "[[", 3)))
      emscan[[9]] <- c(as.numeric(unlist(lapply(dat[4],
                                                "[[", 9))), as.numeric(unlist(lapply(dat[5:length(dat)],
                                                                                     "[[", 4))))
      emscan[[11]] <- c(as.numeric(unlist(lapply(dat[4],
                                                 "[[", 11))), as.numeric(unlist(lapply(dat[5:length(dat)],
                                                                                       "[[", 5))))
      emscan[[14]] <- c(as.numeric(unlist(lapply(dat[4],
                                                 "[[", 14))), as.numeric(unlist(lapply(dat[5:length(dat)],
                                                                                       "[[", 6))))
      emscan
    })
  }
  if(length(emscans) == 1){
    emscans <- unlist(emscans, recursive = FALSE)
  }
  return(emscans)
}
