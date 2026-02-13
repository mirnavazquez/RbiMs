#' Plot pathway-level directional bias (proportion + CI)
#'
#' Plots `prop_pos` with confidence intervals as dot + error bars, and optionally
#' annotates FDR-adjusted p-values.
#'
#' @param pathway_tbl Output tibble from `calc_pathway_directional_bias()`.
#' @param reorder Logical; reorder pathways by `prop_pos`. Default TRUE.
#' @param show_fdr_label Logical; annotate points with FDR text. Default TRUE.
#' @param fdr_digits Digits for scientific formatting. Default 2.
#' @param null_line Numeric; draw a dashed horizontal reference line (e.g. 0.5).
#' @param title Plot title.
#' @param xlab X axis label.
#' @param ylab Y axis label.
#'
#' @return A ggplot object.
#' @export
plot_pathway_directional_bias <- function(pathway_tbl,
                                          reorder = TRUE,
                                          show_fdr_label = TRUE,
                                          fdr_digits = 2,
                                          null_line = 0.5,
                                          title = "Depth bias by pathway (proportion + 95% CI)",
                                          xlab = "rbims_pathway",
                                          ylab = "Proportion of KOs with effect > 0") {
  
  required <- c("rbims_pathway", "prop_pos", "ci_low", "ci_high", "p_adj_fdr")
  missing <- setdiff(required, names(pathway_tbl))
  if (length(missing) > 0) {
    stop("Missing required columns in pathway_tbl: ", paste(missing, collapse = ", "))
  }
  
  df <- pathway_tbl
  
  if (isTRUE(reorder)) {
    df <- df %>%
      dplyr::mutate(rbims_pathway = forcats::fct_reorder(.data$rbims_pathway, .data$prop_pos))
  }
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$rbims_pathway, y = .data$prop_pos)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$ci_low, ymax = .data$ci_high), width = 0.15) +
    ggplot2::geom_hline(yintercept = null_line, linetype = "dashed") +
    ggplot2::theme_bw() +
    ggplot2::labs(title = title, x = xlab, y = ylab)
  
  if (isTRUE(show_fdr_label)) {
    p <- p +
      ggplot2::geom_text(
        ggplot2::aes(
          label = paste0("FDR=", format(.data$p_adj_fdr, scientific = TRUE, digits = fdr_digits))
        ),
        vjust = -1.2,
        size = 3
      )
  }
  
  p
}
