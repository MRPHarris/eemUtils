#' Calculate and extract PARAFAC model residuals
#'
#' @description An expansion upon staRdom::eempf_residuals. Minor internal code changes to fix a couple of bugs, and extra options added.
#'
#' @param pfmodel a PARAFAC model object from staRdom::eem_parafaC()
#' @param eem_list a list of eem objects compliant with the staRdom/eemR framework
#' @param select a character vector containing the names of the desired samples. See ?staRdom::eem_parafac
#' @param cores number of cores used for multi-threading. See ?staRdom::eem_parafac
#' @param denormalise is the PARAFAC model normalised? If so, use the eemlist to denormalise the loadings. Optional; only effects relative magnitude of residuals
#' @param extend_eemlist TRUE/FALSE to ensure equal EEM sizes before calculations. Uses staRdom::eem_extend2largest if differences are found.
#' @param verbose return various messages during operation.
#'
#' @importFrom eemR eem_extract
#' @importFrom tibble column_to_rownames
#' @importFrom dplyr mutate_all
#' @importFrom eemR eem_names
#' @importFrom staRdom eem_extend2largest
#' @importFrom staRdom A_missing
#' @importFrom dplyr mutate
#' @importFrom tidyr gather
#' @importFrom dplyr bind_rows
#'
#' @export
#'
extrpf_residuals <- function(pfmodel, eem_list, select = NULL, cores = parallel::detectCores(logical = FALSE)-1,
                             denormalise = FALSE, extend_eemlist = TRUE, verbose = FALSE, force_names = TRUE){
  cnames <- colnames(pfmodel$A)
  # pfmodel <- norm2A(pfmodel)
  if (!is.null(select)) {
    eem_list <- eem_extract(eem_list, sample = select, keep = TRUE,
                            verbose = FALSE)
  }
  if(isTRUE(extend_eemlist)){
    message("Extending eemlist, if applicable")
    eem_list <- eem_extend2largest(eem_list, interpolation = TRUE)
  }
  # revised normalisation
  if(isTRUE(denormalise)){
    pfmodel$A <- extrpf_loadings_denorm(pfmodel, eemlist = eem_list) %>%
      column_to_rownames(var = 'sample') %>%
      mutate_all(as.numeric)
  }
  if (!all(eem_names(eem_list) %in% rownames(pfmodel$A)) |
      length(eem_list) == 0) {
    if(isTRUE(force_names)){
      rownames(pfmodel$A) <- eem_names(eem_list)
    } else {
      pfmodel <- A_missing(eem_list, pfmodel, cores = cores)
    }
  }
  what <- which(rownames(pfmodel$A) %in% (eem_list %>% eem_names()))
  rnames <- rownames(pfmodel$A)[what]
  pfmodel$A <- as.matrix((pfmodel$A)[what, ]) %>% 'colnames<-'(c(cnames)) %>% 'rownames<-'(c(rnames))
  # building the residuals data. Multiple lapply layers.
  res_data <- lapply(pfmodel$A %>% rownames(), function(sample) {
    # Lapply over samples
    comps <- lapply(pfmodel$A %>% colnames(), function(component) {
      # Another lapply! Prepare component matrices for this sample.
      s1 <- pfmodel$B[,component] %*% t(pfmodel$C[,component])
      s2 <- s1 * pfmodel$A[sample,component]
      s2
    })
    names(comps) <- pfmodel$C %>% colnames() # match excitation band names
    # Combine components to produce fit.
    fit <- comps %>% Reduce("+", .)
    eem <- eem_list[[which(eem_list %>% eem_names == sample)]]
    samp <- eem$x[eem$em %in% rownames(pfmodel$B), eem$ex %in%
                    rownames(pfmodel$C)]
    if(isTRUE(verbose)){
      message(sample)
    }
    res <- samp - fit
    comps <- lapply(pfmodel$A %>% colnames(), function(component) {
      comps[[component]] %>%
        data.frame() %>%
        mutate(type = component,em = rownames(pfmodel$B)) %>%
        gather(ex, value, -em, -type) %>% mutate(ex = substr(ex, 2, 4))
    }) %>% bind_rows()
    colnames(samp) <- rownames(pfmodel$C)
    samp <- samp %>% data.frame() %>% mutate(type = "sample",
                                             em = rownames(pfmodel$B)) %>% gather(ex, value, -em,
                                                                                  -type) %>% mutate(ex = substr(ex, 2, 4))
    rownames(res) <- rownames(pfmodel$B)
    colnames(res) <- rownames(pfmodel$C)
    res <- res %>% data.frame() %>% mutate(type = "residual",
                                           em = rownames(pfmodel$B)) %>% gather(ex, value, -em,
                                                                                -type) %>% mutate(ex = substr(ex, 2, 4))
    res2 <- bind_rows(list(comps, samp, res)) %>% mutate(sample = sample)
  }) %>% bind_rows()
}
