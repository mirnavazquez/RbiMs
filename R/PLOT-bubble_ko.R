#' @title Bubble plot of KO/Pathways/Modules and their values per bin
#' @description Creates a bubble plot of KO/Pathways/Modules per bin.
#' Supports "Abundance", "Binary", "Percentage", or a precomputed table ("None").
#' You can pass extra ggplot2 layers via `...` (e.g., labs(), theme(), guides()).
#'
#' @usage bubble_ko(
#'   tibble_ko, x_axis, y_axis, calc = NULL, data_experiment = NULL,
#'   color_character = NULL, order_bins = NULL, order_metabolism = NULL,
#'   color_pallet = NULL, range_size = NULL, x_labs = TRUE, y_labs = TRUE,
#'   text_x = NULL, text_y = NULL, ...
#' )
#'
#' @param tibble_ko tibble created by mapping_ko/get_subset_* or import functions.
#' @param x_axis bare column name (metabolism table) for the X axis (typically Bin_name).
#' @param y_axis bare column name (metabolism table) for the Y axis (KO/pathway/module).
#' @param calc one of "Abundance", "Binary", "Percentage", or "None".
#' @param data_experiment optional data.frame with metadata (joined by "Bin_name").
#' @param color_character bare column name in metadata/tibble_ko for point color.
#' @param order_bins optional character vector with desired bin order.
#' @param order_metabolism optional character vector with desired metabolism order.
#' @param color_pallet optional character vector of colors for the color scale.
#' @param range_size numeric length-2 vector for point size range (default c(1,5)).
#' @param x_labs logical; if TRUE uses x column name as x label, else NULL.
#' @param y_labs logical; if TRUE uses y column name as y label, else NULL.
#' @param text_x numeric size for x text; default 7.
#' @param text_y numeric size for y text; default 7.
#' @param ... extra ggplot2 layers (e.g., labs(y="hola"), theme(...), guides(...)).
#'
#' @import ggplot2 dplyr tidyr rlang pals
#' @examples
#' \dontrun{
#' bubble_ko(ko_bin_mapp, Bin_name, Module, calc = "Binary",
#'           data_experiment = metadata, color_character = Clades,
#'           labs(y = "hello"), theme(legend.position = "bottom"))
#' }
#' @noRd
bubble_ko <- function(tibble_ko,
                      x_axis,
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
                      text_y = NULL,
                      ...) {
  
  # --- Quosures & labels -------------------------------------------------
  x_axis_enquo <- enquo(x_axis)
  y_axis_enquo <- enquo(y_axis)
  color_character_enquo <- enquo(color_character)
  
  x_axis_label <- as_label(x_axis_enquo)
  y_axis_label <- as_label(y_axis_enquo)
  
  # If X is not Bin_name, swap axes so that X becomes Bin_name branch
  if (x_axis_label != "Bin_name") {
    y_axis_enquo <- enquo(x_axis)
    x_axis_enquo <- enquo(y_axis)
    x_axis_label <- as_label(x_axis_enquo)
    y_axis_label <- as_label(y_axis_enquo)
  }
  
  # --- Defaults ----------------------------------------------------------
  if (is.null(color_pallet)) {
    color_pallet <- as.vector(pals::cols25(20))
    warning("Default color palette of length 20. Provide a longer palette if needed.")
  }
  if (is.null(range_size)) range_size <- c(1, 5)
  if (is.null(text_x)) text_x <- 7
  if (is.null(text_y)) text_y <- 7
  
  # Labels for labs()
  x_lab_val <- if (isTRUE(x_labs)) x_axis_label else NULL
  y_lab_val <- if (isTRUE(y_labs)) y_axis_label else NULL
  
  # --- Compute the plotting metric column `tmp` --------------------------
  if (identical(calc, "Abundance")) {
    tibble_ko_mod <- calc_binary(tibble_ko, !!y_axis_enquo, binary = FALSE) |>
      dplyr::rename(tmp = .data$Abundance)
  } else if (identical(calc, "Binary")) {
    tibble_ko_mod <- calc_binary(tibble_ko, !!y_axis_enquo) |>
      dplyr::rename(tmp = .data$Presence_absence)
  } else if (identical(calc, "Percentage")) {
    tibble_ko_mod <- calc_percentage(tibble_ko, !!y_axis_enquo) |>
      dplyr::rename(tmp = .data$Percentage)
  } else { # "None" -> expect a column to rename
    if ("Presence_absence" %in% names(tibble_ko)) {
      tibble_ko_mod <- dplyr::rename(tibble_ko, tmp = .data$Presence_absence)
    } else if ("Abundance" %in% names(tibble_ko)) {
      tibble_ko_mod <- dplyr::rename(tibble_ko, tmp = .data$Abundance)
    } else if ("Percentage" %in% names(tibble_ko)) {
      tibble_ko_mod <- dplyr::rename(tibble_ko, tmp = .data$Percentage)
    } else {
      stop("Could not find a metric column: need one of 'Presence_absence', 'Abundance', or 'Percentage'.")
    }
  }
  
  Table_with_percentage <- tibble_ko_mod |>
    dplyr::select({{ y_axis_enquo }}, .data$Bin_name, .data$tmp) |>
    tidyr::drop_na() |>
    dplyr::distinct() |>
    dplyr::mutate(
      tmp = as.integer(.data$tmp),
      tmp = dplyr::if_else(.data$tmp == 0L, NA_integer_, .data$tmp)
    )
  
  if (!is.null(data_experiment)) {
    Table_with_percentage <- dplyr::left_join(
      Table_with_percentage, data_experiment, by = "Bin_name"
    )
  }
  
  if (is.null(order_metabolism)) {
    order_metabolism <- Table_with_percentage |>
      dplyr::ungroup() |>
      dplyr::distinct({{ y_axis_enquo }}) |>
      dplyr::pull()
  }
  if (is.null(order_bins)) {
    order_bins <- sort(unique(Table_with_percentage$Bin_name))
  }
  
  # --- Build plot (X is Bin_name branch after potential swap) ------------
  p <- ggplot2::ggplot(
    Table_with_percentage,
    ggplot2::aes(
      x = factor(!!x_axis_enquo, levels = !!order_bins),
      y = factor(!!y_axis_enquo, levels = !!order_metabolism),
      size = .data$tmp,
      color = !!color_character_enquo
    )
  ) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::scale_size(range = range_size)
  
  # Color scale: only apply manual if a color column was provided
  if (!quo_is_null(color_character_enquo)) {
    p <- p + ggplot2::scale_color_manual(values = color_pallet)
  }
  
  p <- p +
    ggplot2::theme_linedraw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = text_x, angle = 45, hjust = 1, vjust = 1),
      axis.text.y = ggplot2::element_text(size = text_y)
    ) +
    ggplot2::labs(x = x_lab_val, y = y_lab_val)
  
  # Legends by calc
  if (identical(calc, "Abundance")) {
    p <- p +
      ggplot2::guides(size = ggplot2::guide_legend(title = "Abundance")) +
      ggplot2::theme(
        legend.title = ggplot2::element_text(size = 16),
        legend.text  = ggplot2::element_text(size = 14)
      )
  } else if (identical(calc, "Binary")) {
    p <- p +
      ggplot2::scale_size_continuous(name = "", labels = "Present") +
      ggplot2::theme(
        legend.title = ggplot2::element_text(size = 16),
        legend.text  = ggplot2::element_text(size = 14)
      )
  } else if (identical(calc, "Percentage")) {
    p <- p +
      ggplot2::guides(size = ggplot2::guide_legend(title = "Percentage")) +
      ggplot2::theme(
        legend.title = ggplot2::element_text(size = 16),
        legend.text  = ggplot2::element_text(size = 14)
      )
  }
  
  # --- Apply extra ggplot2 layers passed in `...` -------------------------
  dots <- rlang::list2(...)
  if (length(dots)) for (lay in dots) p <- p + lay
  
  return(p)
}
