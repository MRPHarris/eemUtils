# These functions relate to the plotting of excitation-emission matrices.

#' Plot an EEM object in 3D using plotly.
#'
#' @description A conversion of staRdom::eempf_comps3D(), applied to EEM objects.
#'
#' @param eem An eem object, compliant with the eemR/staRdom packages.
#' @param origin How was this eem imported? One of 'csv' (generic) or 'DAT' (Aqualog import)
#' @param colour One of "default", "rainbow", or a user-supplied vector of hexadecimal colours.
#'
#' @export
#'
plot_eem_3D <- function (eem, origin = "DAT", colour = "default"){
  if(!(origin == "DAT" || origin == "csv")){
    warning("please specify an origin for the eem: either 'csv' or 'DAT'")
  }
  if(class(eem) == "eem"){
    eem2df_nested <- function (eem, gather = TRUE){
      eem_df <- eem[["x"]] %>% data.frame()
      colnames(eem_df) <- eem[["ex"]]
      rownames(eem_df) <- eem[["em"]]
      if (gather == TRUE){
        eem_df <- eem_df %>% tibble::rownames_to_column("em") %>%
          gather(ex, value, -em)}
      eem_df
    }
    data <- eem %>% eem2df_nested()
  } else(warning("please use an object of class 'eem'"))
  z <- data %>% data.frame() %>% spread(em, value) %>% remove_rownames() %>% column_to_rownames("ex") %>% as.matrix()
  ex <- data$ex %>% unique() %>% as.numeric() %>% as.vector()
  if(origin == "csv"){
    ex <- rev(ex)
  }
  em <- data$em %>% unique() %>% as.numeric() %>% as.vector()
  scene <- list(
    xaxis = list(
      title = "em"),
    yaxis = list(
      title = "ex",
      autorange = "reversed"))
  if(colour == "default"){
    data("eem_palette_12")
    colpal <- eem_palette_12
    if(sum(eem$x < 0, na.rm = TRUE) == 0){
      newpal <- colorRampPalette(c(colpal[2:length(colpal)]))
      colpal <- newpal(12)
    }
  } else if(colour == "rainbow"){
    colpal <- rainbow(12)[9:1]
  } else {
    colpal <- colour
  }
  plotly::plot_ly(x = em, y = ex, z = z, colors = colpal) %>%
    plotly::layout(scene = scene) %>%
    plotly::add_surface()
}

#' Save all the EEMs within an eemlist object as .png images
#'
#' @description loop through an eemlist, saving EEM data as individual .png images.
#'
#' @param eemlist the target eemlist object compliant with the staRdom/eemR framework.
#' @param contour TRUE/FALSE to add contour to EEM.
#' @param output_dir full path to the destination folder. If NULL, images will be sent to the working directory.
#'
#' @export
#'
save_eemlist_pngs <- function(eemlist, contour = FALSE, output_dir = NULL){
  cyn = contour
  if(is.null(output_dir)){
    output_dir <- paste0(getwd(),"/")
  }
  for (i in seq_along(eemlist)){
    EEM_Name = eemlist[[i]][["sample"]]
    print(EEM_Name)
    message(paste0("EEM ",i,"/",length(eemlist)))
    Target_eemlist <- vector("list", 1)
    Target_eemlist[[1]] <- eemlist[[i]]
    png(paste0(output_dir,EEM_Name,".png"), units="cm", width=18, height=14, res=300)
    print(eem_overview_plot(Target_eemlist, contour = cyn)) # change here to output contours in .png EEMs
    dev.off()
  }
}

#' Extract a specific EEM from a set of eemlists
#'
#' @description when given a list of eemlists representing different stages in a processing
#'    framework, extract a single EEM from each eemlist in order to show the processing steps.
#'
#' @param list_of_eemlist a list vector, with each element comprising an eemlist compliant with the staRdom/eemR framework
#' @param which_eem integer - numerically index which eem within the eemlists?
#' @param output_dir NULL or full path to destination directory in which to save a .png image of the EEMs together.
#' @param contour TRUE/FALSE include contours in exported plots, assuming output_dir != NULL
#'
#' @export
#'
extract_procstep_eems <- function(list_of_eemlists, which_eem = 1, output_dir = NULL, contour = TRUE){
  group_eemlist <- vector(mode = "list", length = 6)
  class(group_eemlist) <- "eemlist"
  group_eemlist[[1]] <- list_of_eemlists[[1]][[which_eem]]
  group_eemlist[[2]] <- list_of_eemlists[[2]][[which_eem]]
  group_eemlist[[3]] <- list_of_eemlists[[3]][[which_eem]]
  group_eemlist[[4]] <- list_of_eemlists[[4]][[which_eem]]
  group_eemlist[[5]] <- list_of_eemlists[[5]][[which_eem]]
  group_eemlist[[6]] <- list_of_eemlists[[6]][[which_eem]]
  name_2 <- list_of_eemlists[[2]][[which_eem]][["sample"]]
  name_3 <- list_of_eemlists[[3]][[which_eem]][["sample"]]
  name_4 <- list_of_eemlists[[4]][[which_eem]][["sample"]]
  name_5 <- list_of_eemlists[[5]][[which_eem]][["sample"]]
  name_6 <- list_of_eemlists[[6]][[which_eem]][["sample"]]
  group_eemlist[[2]][["sample"]] <- paste0(name_2,"_masked")
  group_eemlist[[3]][["sample"]] <- paste0(name_3,"_interp")
  group_eemlist[[4]][["sample"]] <- paste0(name_4,"_desct")
  group_eemlist[[5]][["sample"]] <- paste0(name_5,"_desct2")
  group_eemlist[[6]][["sample"]] <- paste0(name_6,"_final")
  if(!(is.null(output_dir))){
    png(paste0(output_dir,group_eemlist[[1]][["sample"]],"processing_steps.png"), units="cm", width=18, height=14, res=300)
    print(eem_overview_plot(group_eemlist, contour = contour)) # change here to output contours in .png EEMs
    dev.off()
  } else(
    message("no output directory defined; png not exported")
  )
  print(eem_overview_plot(group_eemlist[1:6], spp = 6, contour = contour))
}

#' A tweaked EEM plotter, building on ggeem() from staRdom(). For more detailed information, refer to ?staRdom::ggeem()
#'
#' @description An update to staRdom's existing EEM plotter, ggeem. Option for binning values,
#'      along with new colours amongst other things.
#'
#' @param eem An eem object compliant with the staRdom/eemR framework
#' @param bin_vals Three routes here. NULL for no bins. "colpal" to match the number of bins to the length of the colour palette. Alternatively, provide a numeric value for the number of bins you would like.
#' @param colpal Provide a colour palette as a vector of colours. Two defaults are available - "12pal" for a 12-colour palette based on MATLAB's jet scheme, or "rainbow" for ggeem's default scheme.
#' @param contour TRUE/FALSE to plot contour
#' @param interpolate TRUE/FALSE to interpolate. Fairly certain this is defunct; refer to ?staRdom::ggeem()
#' @param redneg refer to ?staRdom::ggeem()
#' @param legend TRUE/FALSE to display legend
#' @param textsize_multiplier a simple numeric multiplier for increasing text size of graphical elements. To help with R's tricky graphics device text scaling when exporting.
#'
#' @import ggplot2 tidyr eemR
#' @importFrom grDevices rainbow
#'
#' @export
#'
ggeem2 <- function(data, fill_max = FALSE, ...) UseMethod("ggeem2")

#' @rdname ggeem2
#' @export
ggeem2.default <- function(data, fill_max=FALSE, ...){
  stop("Data is not of a suitable format!")
}

#' @rdname ggeem2
#' @export
ggeem2.eemlist <- function(data, fill_max = FALSE, eemlist_order = TRUE, ...)
{
  table <- data %>% lapply(as.data.frame) %>% bind_rows()
  if(isTRUE(eemlist_order)){
    table$sample <- table$sample %>%
      factor(levels = table$sample %>% unique())
  }
  ggeem2(table, fill_max = fill_max ,...)
}

#' @rdname ggeem2
#' @export
ggeem2.eem <- function(data, fill_max = FALSE, ...)
{
  table <- data %>% as.data.frame()
  ggeem2(table, fill_max = fill_max, ...)
}

#' @rdname ggeem2
#' @export
ggeem2.parafac <- function(data, fill_max = FALSE, ...)
{
  table <- data %>% eempf_comp_mat() #eem_list
  table <- lapply(table %>% names(),function(name){
    table[[name]] %>% mutate(sample = name)
  }) %>% bind_rows() %>%
    mutate(sample = factor(sample, levels = colnames(data$A)))
  #filename <- paste0('EEM_PARAFAC_components_',suffix,format(Sys.time(), "%Y%m%d_%H%M%S"))
  ggeem2(table,fill_max=fill_max,...)
}

#' @rdname ggeem2
#' @export
ggeem2.data.frame <- function(data,
                              n_eems = 1,
                              fill_max = FALSE,
                              title_text = NULL,
                              bin_vals = NULL,
                              colpal = "12pal",
                              contour = TRUE,
                              interpolate = FALSE,
                              redneg = NULL,
                              legend = TRUE,
                              textsize_multiplier = 1,...){
  # Colpal handling
  if(colpal[1] == "12pal"){
    #colpal <- (function(...)get(data(...,envir = new.env())))(eem_palette_12) # thanks henfiber https://stackoverflow.com/questions/30951204/load-dataset-from-r-package-using-data-assign-it-directly-to-a-variable
    data("eem_palette_12")
    colpal <- eem_palette_12
    message("Using default colour palette")
    if(sum(data$x < 0, na.rm = TRUE) == 0){
      newpal <- colorRampPalette(c(colpal[2:length(colpal)]))
      colpal <- newpal(12)
    }
  } else if(colpal[1] == "rainbow"){
    colpal <- rainbow(75)[53:1]
    message("Using rainbow colour palette")
  } else if(!is.vector(colpal)) {
    stop("Please provide a vector of colours, or use the defaults!")
  }
  if(n_eems == 0){
    stop("n_eems set to 0. Please provide a number of eems")
  } else if(n_eems == 1){
    # one EEM.
    plot <- ggeem2_single(data,
                          n_eems = n_eems,
                          fill_max = fill_max,
                          title_text = title_text,
                          bin_vals = bin_vals,
                          colpal = colpal,
                          contour = contour,
                          interpolate = interpolate,
                          redneg = redneg,
                          legend = legend,
                          textsize_multiplier = textsize_multiplier,...)
    plot
  } else {
    #data <- table
    plot <- ggeem2_multi(data,
                         n_eems = n_eems,
                         fill_max = fill_max,
                         title_text = title_text,
                         bin_vals = bin_vals,
                         colpal = colpal,
                         contour = contour,
                         interpolate = interpolate,
                         redneg = redneg,
                         legend = legend,
                         textsize_multiplier = textsize_multiplier,...)
    plot
  }
}

#' Single plot handler for ggeem2.
#'
#' @description Plot a single EEM for ggeem2. Params inherit from ggeem2.
#'
#' @noRd
#'
ggeem2_single <- function(data,
                          n_eems = n_eems,
                          fill_max = fill_max,
                          title_text = title_text,
                          bin_vals = bin_vals,
                          colpal = colpal,
                          contour = contour,
                          interpolate = interpolate,
                          redneg = redneg,
                          legend = legend,
                          textsize_multiplier = textsize_multiplier,...){
  eem <- data
  eem_constructed <- eemdf_to_eem(eemdf = eem,
                                  file = NULL,
                                  sample = "constructed ggeem2 eem",
                                  location = NULL,
                                  gathered = TRUE)

  if(colpal[1] == "12pal"){
    #colpal <- (function(...)get(data(...,envir = new.env())))(eem_palette_12) # thanks henfiber https://stackoverflow.com/questions/30951204/load-dataset-from-r-package-using-data-assign-it-directly-to-a-variable
    data("eem_palette_12")
    colpal <- eem_palette_12
    message("Using default colour palette")
    if(sum(eem$x < 0, na.rm = TRUE) == 0){
      newpal <- colorRampPalette(c(colpal[2:length(colpal)]))
      colpal <- newpal(12)
    }
  } else if(colpal[1] == "rainbow"){
    colpal <- rainbow(75)[53:1]
    message("Using rainbow colour palette")
  } else if(!is.vector(colpal)) {
    stop("Please provide a vector of colours, or use the defaults!")
  }
  # slit dimensions
  x_slit_min = eem_constructed$ex[2] - eem_constructed$ex[1]
  x_slit_max = eem_constructed$ex[length(eem_constructed$ex)] - eem_constructed$ex[length(eem_constructed$ex)-1]
  y_slit_min = eem_constructed$em[2] - eem_constructed$em[1]
  y_slit_max = eem_constructed$em[length(eem_constructed$em)] - eem_constructed$em[length(eem_constructed$em)-1]
  # panel border rectangle
  rect <- data.frame(
    x = c(min(eem_constructed$ex),min(eem_constructed$ex),max(eem_constructed$ex),max(eem_constructed$ex)),
    y = c(min(eem_constructed$em),max(eem_constructed$em),min(eem_constructed$em),max(eem_constructed$em))
  )
  if(!is.null(redneg)){
    warning("redneg is deprecated and will be ignored! Please use the argument 'colpal = c(rainbow(75)[58],rainbow(75)[51:1])' to produce similar behaviour.")
  }
  if(!is.null(bin_vals)){
    if(bin_vals == "colpal"){
      message("binning vals based on a max EEM intensity of ",max(eem_constructed$x, na.rm = TRUE), " and ",length(colpal)," bins.")
      eem_df <- eemUtils::eem_bin(eem = eem_constructed,
                                  nbins = length(colpal))
    } else {
      if(!is(bin_vals, "numeric")){
        stop("If setting number of bins, please provide a numeric value")
      }
      message("binning vals based on a max EEM intensity of ",max(eem_constructed$x, na.rm = TRUE), " and ",length(colpal)," bins.")
      eem_df <- eemUtils::eem_bin(eem = eem_constructed,
                                  nbins = bin_vals)
    }
  } else {
    eem_df <- eem
  }
  eem_df$value <- as.numeric(as.character(eem_df$value)) # Just in case there are factors carrying over from eem_bin
  eem_df$ex <- as.numeric(eem_df$ex)
  eem_df$em <- as.numeric(eem_df$em)
  table <- as.data.frame(eem_df)
  if(!is.numeric(fill_max)){
    fill_max <- table$value %>% max(na.rm=TRUE)
  }
  # These lines will fail if the raster package is loaded. Fixed by specifying dplyr for select.
  diffs <- table %>%
    dplyr::select(-value) %>%
    gather("spec","wl", -sample) %>%
    group_by(sample,spec) %>%
    unique() %>%
    #arrange(sample,spec,wl) %>%
    #mutate(diffs = wl - lag(wl))
    summarise(slits = diff(wl) %>% n_distinct()) %>%
    .$slits != 1
  # Plotting
  plot <- table %>%
    ggplot(aes(x = ex, y = em, z = value))+
    labs(x = "Excitation (nm)", y = "Emission (nm)")
  if(any(diffs)){
    plot <- plot +
      layer(mapping = aes(colour = value, fill = value),
            geom = "tile", stat = "identity", position = "identity")
  } else {
    plot <- plot +
      layer(mapping = aes(fill = value),
            geom = "raster", stat = "identity", position = "identity")
  }
  plot <- plot +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    facet_wrap(~ sample)
  if(contour){
    plot <- plot +
      geom_contour(colour = "white", size = 0.2, alpha = 0.5)
  }
  # Binned? Discrete scale.
  if(isTRUE(bin_vals)){
    if(table$value %>% min(na.rm=TRUE) < 0){
      ## Vals for breaks if using bins. Lots of lines because this was a nightmare to wrap my tiny brain around.
      vals_test <- unique(table$value, na.rm = TRUE)[which(!is.na(unique(table$value, na.rm = TRUE)))]
      vals_test_break = append(diff(vals_test),diff(vals_test)[length(diff(vals_test))], after = length(diff(vals_test)))
      vals_shift <- vals_test - (vals_test_break/2)
      vals_shift <- vals_shift[-1]
      vals_test_endbreak_half <- (vals_test_break/2)[length(vals_test_break)]
      vals_shift <- append(vals_shift, (vals_shift[length(vals_shift)] + vals_test_endbreak_half), after = length(vals_shift))
      vals_shift <- round(vals_shift,2)
      vals_labels <- round(vals_test,2)
      # Adaptive shifting in case of length mismatch.
      if(length(colpal) < length(vals_labels)){
        vals_labels <- vals_labels[-1]
        vals_shift <- vals_shift[-length(vals_shift)]
      }
      plot <- plot +
        scale_fill_stepsn(colours = colpal, breaks = vals_shift, labels = vals_labels, limits = c(table$value %>% min(na.rm=TRUE),fill_max),
                          na.value="white") +
        scale_colour_stepsn(colours = colpal, breaks = vals_shift,labels = vals_labels, limits = c(table$value %>% min(na.rm=TRUE),fill_max),
                            na.value="white") +
        scale_x_continuous(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0))
    } else {
      plot <- plot +
        scale_fill_stepsn(colours = colpal, limits = c(0,fill_max), na.value="white")+
        scale_colour_stepsn(colours = colpal, limits = c(0,fill_max), na.value="white") +
        scale_x_continuous(expand = c(0,0))+#, limits = c(min(eem$ex))) +
        scale_y_continuous(expand = c(0,0))#, limits = c(min(eem$em)))
    }
  } else {
    # Not binned - continuous scale.
    if(table$value %>% min(na.rm=TRUE) < 0){
      vals <- c(table$value %>% min(na.rm = TRUE), seq(from = 0, to = fill_max, length.out = length(colpal) - 1))
      vals <- (vals - min(vals))/diff(range(vals))
      plot <- plot +
        scale_fill_gradientn(colours = colpal, values = vals, limits = c(table$value %>% min(na.rm=TRUE),fill_max),
                             na.value="white")+
        scale_colour_gradientn(colours = colpal, values = vals, limits = c(table$value %>% min(na.rm=TRUE),fill_max),
                               na.value="white") +
        scale_x_continuous(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0))
    } else {
      plot <- plot +
        scale_fill_gradientn(colours = colpal, limits = c(0,fill_max), na.value="white")+
        scale_colour_gradientn(colours = colpal, limits = c(0,fill_max), na.value="white") +
        scale_x_continuous(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0))
    }
  }
  # Title handling.
  if(!is.null(title_text)){
    plot <- plot +
      labs(title = title_text) +
      theme(
        plot.title = element_text(hjust = 0.5,
                                  size = 11*textsize_multiplier,
                                  face = "bold")
      )
  } else {
    plot <- plot +
      theme(
        plot.title = element_blank()
      )
  }
  # Final thematic elements.
  plot <- plot +
    theme(
      panel.background = element_rect(fill = 'white', colour = 'black'),
      strip.background = element_blank(),
      strip.text = element_blank(),
      axis.text = element_text(size = 9*textsize_multiplier),
      axis.title = element_text(size = 10*textsize_multiplier),
      legend.text = element_text(size = 9*textsize_multiplier),
      legend.title = element_text(size = 10*textsize_multiplier)
    )
  # Border rectangle
  plot  <- plot +
    annotate(geom = "rect", xmin = min(eem_constructed$ex)-(x_slit_min/2), xmax = max(eem_constructed$ex)+(x_slit_max/2), ymin = min(eem_constructed$em)-(y_slit_min/2), ymax = max(eem_constructed$em)+(y_slit_max/2),
             colour = "black", fill = "white", alpha = 0, size = 0.5)
  # Legend removal
  if(!isTRUE(legend)){
    plot <- plot +
      theme(
        legend.position = "none"
      )
  } #else {   # Unhash this code if you would like a nested legend.
  #plot <- plot +
  #  theme(
  #    legend.position = c(0.9,0.2)
  #  )
  #}
  plot

}

#' Multi-plot handler for ggeem2.
#'
#' @description Plot a single EEM for ggeem2. Params inherit from ggeem2.
#'
#' @noRd
#'
ggeem2_multi <- function(data,
                         n_eems = n_eems,
                         fill_max = fill_max,
                         title_text = title_text,
                         bin_vals = bin_vals,
                         colpal = colpal,
                         contour = contour,
                         interpolate = interpolate,
                         redneg = redneg,
                         legend = legend,
                         textsize_multiplier = textsize_multiplier,...){
  if(!is.null(bin_vals)){
    message("Binning not supported for multiple EEMs with ggeem2. Skipping binning.")
  }
  if(!is.null(redneg)){
    warning("redneg is deprecated and will be ignored! Please use the argument 'colpal = c(rainbow(75)[58],rainbow(75)[51:1])' to produce similar behaviour.")
  }
  table <- data %>%
    mutate_at(vars(ex,em,value),as.numeric)
  if(!is.numeric(fill_max)){
    fill_max <- table$value %>% max(na.rm=TRUE)
  }
  diffs <- table %>%
    select(-value) %>%
    gather("spec","wl", -sample) %>%
    group_by(sample,spec) %>%
    unique() %>%
    summarise(slits = diff(wl) %>% n_distinct()) %>%
    .$slits != 1

  plot <- table %>%
    ggplot(aes(x = ex, y = em, z = value))+
    labs(x = "Excitation (nm)", y = "Emission (nm)")

  if(any(diffs)){
    plot <- plot +
      layer(mapping = aes(colour = value, fill = value),
            geom = "tile", stat = "identity", position = "identity")
  } else {
    plot <- plot +
      layer(mapping = aes(fill = value),
            geom = "raster", stat = "identity", position = "identity")
  }

  plot <- plot +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    facet_wrap(~ sample)
  if(contour){
    plot <- plot +
      geom_contour(colour = "black", size = 0.3)
  }
  if(table$value %>% min(na.rm=TRUE) < 0){
    vals <- c(table$value %>% min(na.rm = TRUE), seq(from = 0, to = fill_max, length.out = length(colpal) - 1))
    vals <- (vals - min(vals))/diff(range(vals))
    plot <- plot +
      scale_fill_gradientn(colours = colpal, values = vals, limits = c(table$value %>% min(na.rm=TRUE),fill_max),
                           na.value="white")+
      scale_colour_gradientn(colours = colpal, values = vals, limits = c(table$value %>% min(na.rm=TRUE),fill_max),
                             na.value="white") +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0))
  } else {
    plot <- plot +
      scale_fill_gradientn(colours = colpal, limits = c(0,fill_max), na.value="white")+
      scale_colour_gradientn(colours = colpal, limits = c(0,fill_max), na.value="white") +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0))
  }
  if(!is.null(title_text)){
    plot <- plot +
      labs(title = title_text) +
      theme(plot.title = element_text(hjust = 0.5, size = 11*textsize_multiplier, face = "bold"))
  } else {
    plot <- plot +
      theme(plot.title = element_blank())
  }
  plot <- plot +
    theme(
      panel.background = element_rect(fill = 'white', colour = 'black'),
      strip.background = element_blank(),
      strip.text = element_blank(),
      axis.text = element_text(size = 9*textsize_multiplier),
      axis.title = element_text(size = 10*textsize_multiplier),
      legend.text = element_text(size = 9*textsize_multiplier),
      legend.title = element_text(size = 10*textsize_multiplier)
    )
  #plot  <- plot +
  #  annotate(geom = "rect", xmin = min(eem_constructed$ex)-(x_slit_min/2), xmax = max(eem_constructed$ex)+(x_slit_max/2), ymin = min(eem_constructed$em)-(y_slit_min/2), ymax = max(eem_constructed$em)+(y_slit_max/2),
  #           colour = "black", fill = "white", alpha = 0, size = 0.5)
  if(!isTRUE(legend)){
    plot <- plot +
      theme(legend.position = "none")
  }
  plot
}
