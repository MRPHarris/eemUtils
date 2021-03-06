% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eem_data_functions.R
\name{eem_setNA_mod}
\alias{eem_setNA_mod}
\title{Set areas of one or more EEMs to NA.}
\usage{
eem_setNA_mod(
  eemlist,
  sample = NULL,
  em = NULL,
  ex = NULL,
  interpolate = FALSE
)
}
\arguments{
\item{eemlist}{A list of EEMs in a format compliant with eemR/staRdom.}

\item{sample}{optional, names or indices of samples to process}

\item{em}{optional, emission wavelengths to set NA}

\item{ex}{optional, excitation wavelengths to set NA}

\item{interpolate}{FALSE, 1 or 2, interpolate NAs or not, 2 different methods.}
}
\description{
A tweaked/re-built version of staRdom::set_NA(). Sets areas of
       eems to NA. Option to interpolate. Changes from set_NA() relate to how
       the range is identified. I recommend interpolating seperately.
}
