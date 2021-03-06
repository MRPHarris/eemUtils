% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eem_data_functions.R
\name{eem_0_to_NA}
\alias{eem_0_to_NA}
\title{Set all instances of '0'/zero to NA within a list of EEMs.}
\usage{
eem_0_to_NA(eemlist, outputfolder = NULL)
}
\arguments{
\item{eemlist}{A list of EEMs in a format compliant with eemR/staRdom.}

\item{outputfolder}{Path to the folder where the new eemlist will be sent. Or NULL, if no export is desired.}
}
\description{
Go through a list of eemR/staRdom-compliant EEMs, and replace all
      instances of 0 with NA. Used in cases where e.g. an Aqualog Rayleigh-masking
      step has masked those rows with 0, preventing proper interpolation and
      PARAFAC modeling.
}
