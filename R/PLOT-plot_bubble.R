#' @title Bubble plot of features (KO/pathways/modules) per bin
#' @description Generic bubble plot for RBIMS tables (KEGG/InterPro/Pfam/...).
#' Supports "Abundance", "Binary", "Percentage", or "None" (precomputed).
#' You can set axis labels via x_labs/y_labs directly (TRUE/FALSE/character).
#' Extra ggplot2 layers can be passed via `...` (theme(), guides(), etc.).
#'
#' @param tibble_ko Tibble produced by mapping/get_subset/import functions.
#' @param x_axis Bare column (usually Bin_name).
#' @param y_axis Bare column (KO/pathway/module/etc.).
#' @param analysis Optional string ("KEGG", "InterPro", etc.). Informative only.
#' @param calc "Abundance", "Binary", "Percentage", or "None".
#' @param data_experiment Optional metadata (joined by "Bin_name").
#' @param color_character Bare column for color (metadata or tibble_ko).
#' @param order_bins Optional character vector with desired bin order.
#' @param order_metabolism Optional character vector with desired y order.
#' @param color_pallet Optional character vector for manual color scale.
#' @param range_size Numeric length-2 for point size range (default c(1,5)).
#' @param x_labs TRUE uses x name, FALSE uses NULL, character sets label.
#' @param y_labs TRUE uses y name, FALSE uses NULL, character sets label.
#' @param text_x Numeric axis text size (default 7).
#' @param text_y Numeric axis text size (default 7).
#' @param ... Extra ggplot2 layers (e.g., theme(...), guides(...)).
#'
#' @import ggplot2 dplyr tidyr rlang pals
#' @export
plot_bubble <- function(tibble_ko,
                        x_axis,
                        y_axis,
                        analysis = NULL,
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
  
  # ---- Early guard to satisfy tests ------------------------------------
  # If calc is missing AND tibble_ko doesn't have a metric column, throw the
  # exact message expected by tests.
  has_metric_col <- any(c("Presence_absence", "Abundance", "Percentage") %in% names(tibble_ko))
  if (is.null(calc) && !has_metric_col) {
    stop("calc must have a value between Abundance, Binary or Percentage")
  }
  
  # ---- Quosures & labels ------------------------------------------------
  x_axis_enquo <- rlang::enquo(x_axis)
  y_axis_enquo <- rlang::enquo(y_axis)
  color_character_enquo <- rlang::enquo(color_character)
  
  x_axis_label <- rlang::as_label(x_axis_enquo)
  y_axis_label <- rlang::as_label(y_axis_enquo)
  
  # If X is not Bin_name, swap to keep plotting branch consistent
  if (x_axis_label != "Bin_name") {
    y_axis_enquo <- rlang::enquo(x_axis)
    x_axis_enquo <- rlang::enquo(y_axis)
    x_axis_label <- rlang::as_label(x_axis_enquo)
    y_axis_label <- rlang::as_label(y_axis_enquo)
  }
  
  # ---- Defaults ---------------------------------------------------------
  if (is.null(range_size)) range_size <- c(1, 5)
  if (is.null(text_x)) text_x <- 7
  if (is.null(text_y)) text_y <- 7
  
  # Labels: allow TRUE/FALSE/character
  lbl <- function(val, fallback) if (isTRUE(val)) fallback else if (is.character(val)) val else NULL
  x_lab_val <- lbl(x_labs, x_axis_label)
  y_lab_val <- lbl(y_labs, y_axis_label)
  
  # ---- Compute metric column `tmp` --------------------------------------
  if (identical(calc, "Abundance")) {
    tibble_ko_mod <- calc_binary(tibble_ko, !!y_axis_enquo, binary = FALSE) |>
      dplyr::rename(tmp = .data$Abundance)
  } else if (identical(calc, "Binary")) {
    tibble_ko_mod <- calc_binary(tibble_ko, !!y_axis_enquo) |>
      dplyr::rename(tmp = .data$Presence_absence)
  } else if (identical(calc, "Percentage")) {
    tibble_ko_mod <- calc_percentage(tibble_ko, !!y_axis_enquo) |>
      dplyr::rename(tmp = .data$Percentage)
  } else { # "None": rename an existing metric column
    if ("Presence_absence" %in% names(tibble_ko)) {
      tibble_ko_mod <- dplyr::rename(tibble_ko, tmp = .data$Presence_absence)
    } else if ("Abundance" %in% names(tibble_ko)) {
      tibble_ko_mod <- dplyr::rename(tibble_ko, tmp = .data$Abundance)
    } else if ("Percentage" %in% names(tibble_ko)) {
      tibble_ko_mod <- dplyr::rename(tibble_ko, tmp = .data$Percentage)
    } else {
      # fall back to the same test message for consistency
      stop("calc must have a value between Abundance, Binary or Percentage")
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
  
  # ---- Join metadata if provided ---------------------------------------
  if (!is.null(data_experiment)) {
    Table_with_percentage <- dplyr::left_join(
      Table_with_percentage, data_experiment, by = "Bin_name"
    )
  }
  
  # ---- Orders (bins & y) ------------------------------------------------
  if (is.null(order_metabolism)) {
    order_metabolism <- Table_with_percentage |>
      dplyr::ungroup() |>
      dplyr::distinct({{ y_axis_enquo }}) |>
      dplyr::pull()
  }
  if (is.null(order_bins)) {
    order_bins <- sort(unique(Table_with_percentage$Bin_name))
  }
  
  # ---- Color palette logic ---------------------------------------------
  # Only set (and possibly warn) if color aesthetic is used
  use_color <- rlang::as_label(color_character_enquo) != "NULL"
  if (use_color && is.null(color_pallet)) {
    # compute number of distinct color groups if possible
    n_cols <- tryCatch({
      n <- Table_with_percentage |>
        dplyr::pull(!!color_character_enquo) |>
        unique() |>
        length()
      n
    }, error = function(e) NA_integer_)
    # default palette
    color_pallet <- as.vector(pals::cols25(20))
    # warn only if we can estimate that palette is too short
    if (!is.na(n_cols) && n_cols > length(color_pallet)) {
      warning("Default color palette (n=20). Provide a larger palette if needed.")
    }
  }
  
  # ---- Build plot -------------------------------------------------------
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
    ggplot2::scale_size(range = range_size) +
    ggplot2::theme_linedraw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = text_x, angle = 45, hjust = 1, vjust = 1),
      axis.text.y = ggplot2::element_text(size = text_y)
    ) +
    ggplot2::labs(x = x_lab_val, y = y_lab_val)
  
  if (use_color) {
    p <- p + ggplot2::scale_color_manual(values = color_pallet)
  }
  
  # ---- Extra ggplot2 layers from `...` ---------------------------------
  dots <- rlang::list2(...)
  if (length(dots)) {
    for (lay in dots) {
      # be conservative: add only common ggplot components (omit 'labels' to avoid <labels> trap)
      if (inherits(lay, c("gg", "ggplot", "LayerInstance", "Scale", "Coord", "Facet", "theme"))) {
        p <- p + lay
      }
      # silently ignore anything else (including 'labels')
    }
  }
  
  # ---- Legends by calc --------------------------------------------------
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
  
  return(p)
}
