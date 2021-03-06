% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc_data_functions.R
\name{ssc_plot}
\alias{ssc_plot}
\title{Plot two spectra and show congruence metric data}
\usage{
ssc_plot(mat1, mat2, normalise = TRUE, labels = NULL)
}
\arguments{
\item{mat1}{a matrix}

\item{mat2}{a matrix}

\item{normalise}{TRUE/FALSE to normalise the matrices before comparison and subsequent plotting.}

\item{labels}{Optional character vectors to label the spectra with. Correspond to mat1, mat2.}
}
\description{
plots two matrices (e.g. spectra) in a similar fashion to figure 2 from Wunsch et al., 2019.
}
