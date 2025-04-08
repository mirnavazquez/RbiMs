#' @title Bubble plot of MEROPS and its relative abundance 
#' within each bin.
#' @description Creates a bubble plot of MEROPS relative 
#' abundance within each bin. 
#' It uses metadata information to color bubbles.
#' @usage bubble_merops(tibble_ko, x_axis, y_axis, calc = NULL,
#' data_experiment=NULL, color_character=NULL, order_bins=NULL, 
#' order_metabolism=NULL, color_pallet=NULL, range_size=NULL, 
#' x_labs=TRUE, y_labs=TRUE, text_x=NULL, text_y=NULL)
#' @param tibble_ko A tibble object created with the mapping_ko 
#' or get_subset_* functions. 
#' @param x_axis A string, column name for the x-axis. 
#' Usually "Bin_name" or a metabolic category.
#' @param y_axis A string, column name for the y-axis. 
#' Usually a MEROPS or PFAM identifier.
#' @param calc Optional. Either "Abundance" or "Binary". Defines which 
#' measure to calculate (relative abundance or presence/absence).
#' @param data_experiment Optional. A data frame with metadata info.
#' @param color_character Optional. A column name from metadata 
#' used to color bubbles.
#' @param order_bins Optional. Character vector defining bin order.
#' @param order_metabolism Optional. Character vector for metabolism order.
#' @param color_pallet Optional. Character vector of color hex codes.
#' @param range_size Optional. Numeric vector for point size range.
#' @param x_labs Optional. TRUE for x-axis labels, FALSE for none.
#' @param y_labs Optional. TRUE for y-axis labels, FALSE for none.
#' @param text_x Optional. Numeric value for x-axis text size.
#' @param text_y Optional. Numeric value for y-axis text size.
#' @details This function is part of a package for bin metabolism analysis.
#' @import ggplot2 dplyr rlang pals
#' @examples
#' \dontrun{
#' bubble_merops(tibble_ko=input_data, x_axis=Bin_name, y_axis=PFAM, 
#' data_experiment=metadata, color_character=Genus)
#' }
#' @noRd
bubble_merops <- function(tibble_ko,
                          x_axis = Bin_name, 
                          y_axis,
                          calc = NULL,
                          data_experiment = NULL,
                          color_character = NULL,
                          order_bins = NULL,
                          order_metabolism = NULL,
                          color_pallet = NULL, 
                          range_size = NULL,
                          x_labs = TRUE,
                          y_labs = TRUE,
                          text_x = NULL,
                          text_y = NULL) {
  
  # Enquote inputs
  x_axis_enquo <- enquo(x_axis)
  y_axis_enquo <- enquo(y_axis)
  x_axis_label <- as_label(x_axis_enquo)
  y_axis_label <- as_label(y_axis_enquo)
  color_character_enquo <- enquo(color_character)
  
  # Flip axes if Bin_name is not the x-axis
  if (x_axis_label != "Bin_name") {
    y_axis_enquo <- enquo(x_axis)
    x_axis_enquo <- enquo(y_axis)
    x_axis_label <- as_label(x_axis_enquo)
    y_axis_label <- as_label(y_axis_enquo)
  }
  
  # Set default color palette
  if (is.null(color_pallet)) {
    color_pallet <- as.vector(cols25(20))
  }
  
  # Set default size range
  if (is.null(range_size)) {
    range_size <- c(1, 10)
  }
  
  # Axis labels
  x_labs <- if (isTRUE(x_labs)) x_axis_label else NULL
  y_labs <- if (isTRUE(y_labs)) y_axis_label else NULL
  
  # Axis text size
  if (is.null(text_x)) text_x <- 7
  if (is.null(text_y)) text_y <- 7
  
  # Calculate abundance or binary data
  if (calc == "Abundance") {
    tibble_ko_mod <- calc_binary(tibble_ko, !!y_axis_enquo, binary = FALSE) %>%
      rename(tmp = .data$Abundance)
  } else if (calc == "Binary") {
    tibble_ko_mod <- calc_binary(tibble_ko, !!y_axis_enquo, binary = TRUE) %>%
      rename(tmp = .data$Presence_absence)
  } else {
    tibble_ko_mod <- tibble_ko
    tibble_ko_mod$tmp <- 1  # fallback: assume presence
  }
  
  # Build bubble data
  bubble <- tibble_ko_mod %>%
    dplyr::select(!!y_axis_enquo, Bin_name, tmp) %>%
    drop_na() %>%
    distinct() %>%
    mutate(tmp = ifelse(tmp == 0, NA_integer_, as.integer(tmp)))
  
  # Set metabolism order
  if (is.null(order_metabolism)) {
    order_metabolism <- bubble %>%
      ungroup() %>%
      dplyr::select(!!y_axis_enquo) %>%
      distinct() %>%
      pull()
  }
  
  # Set bin order
  if (is.null(order_bins)) {
    order_bins <- sort(unique(bubble$Bin_name))
  }
  
  # Re-enquote for plotting
  x_axis_enquo <- enquo(x_axis)
  y_axis_enquo <- enquo(y_axis)
  
  # Create plot
  if (x_axis_label == "Bin_name") {
    plot_bubble <- ggplot(bubble, aes(
      x = factor(!!x_axis_enquo, levels = order_bins),
      y = factor(!!y_axis_enquo, levels = order_metabolism),
      size = tmp,
      color = !!color_character_enquo)) +
      geom_point(alpha = 0.5) +
      scale_size(range = range_size) +
      scale_color_manual(values = color_pallet) +
      theme_linedraw() +
      theme(
        axis.text.x = element_text(size = text_x, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = text_y)) +
      xlab(x_labs) +
      ylab(y_labs)
  } else {
    plot_bubble <- ggplot(bubble, aes(
      x = factor(!!x_axis_enquo, levels = order_metabolism),
      y = factor(!!y_axis_enquo, levels = order_bins),
      size = tmp,
      color = !!color_character_enquo)) +
      geom_point(alpha = 0.5) +
      scale_size(range = range_size) +
      scale_color_manual(values = color_pallet) +
      theme_linedraw() +
      theme(
        axis.text.x = element_text(size = text_x, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = text_y)) +
      xlab(x_labs) +
      ylab(y_labs)
  }
  
  # Final plot adjustments by calc type
  if (calc == "Abundance") {
    bubble <- bubble %>% rename(Abundance = tmp)
    plot_bubble <- plot_bubble + guides(size = guide_legend(title = "Abundance"))
  } else if (calc == "Binary") {
    bubble <- bubble %>% rename("Presence/Absence" = tmp)
    plot_bubble <- plot_bubble +
      scale_size_continuous(name = "", labels = "Present", range = range_size) +
      theme(
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)
      )
  }
  
  suppressWarnings(return(plot_bubble))
}
