#' Get SSC along with alpha and beta penalty terms
#'
#' @description The shift- and shape sensitive congruence (SSC) was developed by Wunsch et al., 2019 as an improvement
#'      upon the TCC metric. It incorporates two penalty terms, alpha and beta, to account for differences in the wavelength
#'      peak position and area. This function adds these terms to the data frame returned by staRdom::ssc()
#'
#' @param mat1 a matrix
#' @param mat2 a matrix
#' @param tcc TRUE/FALSE to return only TCC value instead of SSC, alpha and beta.
#'
#' @export
#'
ssc_more <- function (mat1, mat2, tcc = FALSE) {
  if (any(is.null(mat1), is.na(mat1), is.null(mat2), is.na(mat2))) {
    a <- NA
  } else {
    a <- lapply(1:ncol(mat1), function(nc) {
      col1 <- mat1[, nc]
      apply(mat2, 2, function(col2) {
        tcc_cal <- sum(col1 * col2)/sqrt(sum(col1^2) *
                                           sum(col2^2))
        if (!tcc) {
          wl <- as.numeric(names(col1))
          if (any(is.na(wl)) | pracma::isempty(wl)) {
            stop("SSCs cannot be calculated. Please add wavelengths as rownames of the matrices!")
          }
          alpha <- abs((wl[which.max(col1)] - wl[which.max(col2)])/diff(range(wl)))
          beta <- abs((sum(col1/max(col1)) - sum(col2/max(col2)))/diff(range(wl)))
          ssc <- tcc_cal - alpha - beta
          ssc <- c(ssc, alpha, beta)
          # rownames(a) <- (c("ssc","alpha","beta"))
        } else {
          tcc_cal
        }
      })
    }) %>% setNames(colnames(mat1)) %>% do.call(rbind, .)
  }
  attr(a, "method") <- ifelse(tcc, "TCC", "SSC")
  if(!isTRUE(tcc)){
    rownames (a) <- c("ssc","alpha","beta")
  }
  a
}
