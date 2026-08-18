#' Calculate the sample loadings for samples not involved in model building
#'
#' @description This is simply a port of staRdom's A_missing function, with an option to skip reduction of EEMs.
#'
#' @param eem_list object of class eemlist with sample data
#' @param pfmodel object of class parafac
#' @param skip_reduce TRUE/FALSE to force reduction of EEM dimensions to their smallest member. Required for modelling, but can cause errors in practice.
#' @param cores number of cores to use for parallel processing
#' @param components optionally supply components to use manually, either as a variable of class parafac_components or as a list of variables of class parafac_components, if you do so,
#' @param const optional constraints for model, just used, when components are supplied
#' @param control optional constraint control parameters for model, just used, when components are supplied
#' @param ... additional arguments passed to eem_parafac
#'
#' @importFrom staRdom eem_exclude
#' @importFrom staRdom eempf_bindxc
#' @importFrom staRdom eem_red2smallest
#' @importFrom staRdom eem_parafac
#'
#' @export
#'
A_missing_mod <- function (eem_list,
                           pfmodel = NULL,
                           skip_reduce = TRUE,
                           cores = parallel::detectCores(logical = FALSE),
                           components = NULL,
                           const = NULL,
                           control = NULL, ...)
{
  # Previous use of eem_red2smallest caused a range of errors; option to skip.
  if(!isTRUE(skip_reduce)){
    eem_list <- eem_red2smallest(eem_list)
  }
  if (is.null(pfmodel) & is.null(components)){
    stop("You must either specify a model or components as a base for the newly generated model!")
  }
  exclude <- list(ex = eem_list[[1]]$ex[!(eem_list[[1]]$ex %in% rownames(pfmodel$C))],
                  em = eem_list[[1]]$em[!(eem_list[[1]]$em %in% rownames(pfmodel$B))],
                  sample = c())
  x <- eem_list %>%
    eem_exclude(exclude)
  if (!is.null(components)) {
    if (!is.null(pfmodel))
      warning("The base model is ignored since you provided components manually!")
    if (inherits(components[[1]], "parafac_components"))
      components <- eempf_bindxc(components)
    if (inherits(components, "parafac_components")) {
      Bfixed <- components$B
      Cfixed <- components$C
      comps <- ncol(components$B)
      if (is.null(const))
        const <- c("nonneg", "nonneg", "nonneg")
    } else {
      stop("The list of components you supplied is invalid!")
    }
  } else {
    Bfixed <- pfmodel$B
    Cfixed <- pfmodel$C
    comps <- pfmodel$A %>% ncol()
    normalise = (!is.null(attr(pfmodel, "norm_factors")))
    control = pfmodel$control
    const = pfmodel$const
  }
  missingAs <- eem_parafac(x, comps = comps, normalise = (!is.null(attr(pfmodel,
                                                                        "norm_factors"))), Bfixed = Bfixed, Cfixed = Cfixed,
                           cores = cores, const = const, control = control, ...)
  missingAs[[1]]
}
