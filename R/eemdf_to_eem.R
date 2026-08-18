#' Convert a dataframe to an EEM object of the style used by EEM/eemR/staRdom.
#'
#' @description An alternative to eemR's eem constructor. Intended to reverse as.data.frame(eem).
#'
#' @param eemdf the dataframe to be coerced to an EEM object.
#' @param file filename of the EEM, if applicable.
#' @param sample the samplename of the EEM, if applicable.
#' @param location the location of the EEM file, if applicable.
#' @param gathered TRUE/FALSE is the eemdf in a short (not gathered; FALSE) or a long (gathered; TRUE) format?
#'
#' @export
#'
eemdf_to_eem <- function(eemdf,
                         file = NULL,
                         sample = NULL,
                         location = NULL,
                         gathered = FALSE){
  # code adapted from staRdom's .eem_csv importer.
  x <- eemdf
  if(!isTRUE(gathered)){
    # The eem is in a short, non-gathered format.
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
  } else {
    # The eem is in a long or 'gathered' format.
    gath_df <- eemdf
    gath_df$value <- as.numeric(as.character(gath_df$value)) # factor management
    if("sample" %in% colnames(gath_df)){
      gath_df <- subset(gath_df, select = -c(sample))
    }
    if(colnames(gath_df)[3] == "z"){
      gath_df_short <- spread(data = gath_df, key = "ex", value = "z")
    } else {
      gath_df_short <- spread(data = gath_df, key = "ex", value = "value")
    }
    rnames <- as.matrix(gath_df_short[,1])
    gath_df_short <- as.data.frame(gath_df_short) %>%
      'rownames<-'(c(rnames)) %>%
      select(-c(1))
    eem <- eemdf_to_eem(eemdf = gath_df_short,
                        file = file,
                        sample = sample,
                        location = location,
                        gathered = FALSE) %>% 'class<-'(c('eem'))
    return(eem)
  }
}
