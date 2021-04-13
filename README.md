
<!-- README.md is generated from README.Rmd. Please edit that file -->

# eemUtils

<!-- badges: start -->
<!-- badges: end -->

The **eemUtils** package serves as a repository for various functions to
manipulate, analyse, and plot fluorescence data. It utilises the
existing framework provided by the packages
[eemR](https://cran.r-project.org/web/packages/eemR/index.html) and
[staRdom](https://github.com/MatthiasPucher/staRdom). Many of the
functions within this package are alterations to existing functions
within from eemR or staRdom - thus, if you use the functions from this
package, please allocate proper credit to those packages and their
authors. Most eemUtils functions deal with fluorescence data, though
some are general utility functions which, whilst a part of the
fluorescence functions, can be used more generally.

This package is a work in progress - many of the functions were
developed for a dataset currently in preparation for publication, and so
I cannot currently upload data or detailed examples. The goal is to
eventually upload that dataset and the accompanying workflow as an
offshoot to this package.

If you have any questions or comments, please don’t hesitate to get in
touch.

## Example functions

Below are just a few of the 20+ functions contained within this package.

**plot\_eem\_3D**

A simple conversion of `staRdom::eempf_comps3D()`, but for use with
sample EEM data, rather than outputs from a PARAFAC model.

    plot_eem_3D()

<p align="center">
<img src="man/figures/3D_eem_example.png" height="400px" />
</p>

**Generate\_CORCONDIA**

`Generate_CORCONDIA()` is a simple function wrapper for staRdom’s
existing core consistency diagnostic function
`staRdom::eempf_corcondia()`. It produces a more legible output.

<p align="center">
<img src="man/figures/generate_CORCONDIA_example.PNG" height="400px" />
</p>

**extrpf\_loadings**

Extract the modeled per-sample fluorescence intensity loadings for each
component within a PARAFAC model. This is a simple way to get quick
series data from any number of PARAFAC models generated by
`staRdom::eem_parafac()`

## Installation

To get access to the functions in **eemUtils**, simply use the
**devtools** package to install the package from github.

``` r
devtools::install_github("MRPHarris/eemUtils")
```
