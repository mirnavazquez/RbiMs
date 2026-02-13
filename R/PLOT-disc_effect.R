#' Plot differential effect for discriminant features
#'
#' This plot shows the effect size (x-axis) for top discriminant features.
#' Point size reflects significance as -log10(FDR), and color encodes effect direction.
#'
#' @param disc_obj Object stored as attr(, "rbims_disc") from get_subset_discriminant().
#' @param top_n Number of features to display.
#'
#' @return A ggplot object.
#' @export
plot_disc_effect <- function(disc_obj, top_n = 40) {
  
  cons <- disc_obj$consensus %>%
    dplyr::filter(!is.na(.data$effect)) %>%
    dplyr::mutate(
      neg_logFDR = -log10(.data$aldex_FDR)
    ) %>%
    dplyr::arrange(dplyr::desc(.data$score),
                   .data$aldex_FDR,
                   dplyr::desc(abs(.data$effect))) %>%
    dplyr::slice_head(n = top_n)
  
  ggplot2::ggplot(
    cons,
    ggplot2::aes(
      x     = .data$effect,
      y     = forcats::fct_reorder(.data$feature, .data$effect),
      size  = .data$neg_logFDR,
      color = .data$effect
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::scale_color_gradient2(
      low  = "dodgerblue3",
      mid  = "grey70",
      high = "firebrick3",
      name = "Effect"
    ) +
    ggplot2::scale_size_continuous(name = "-log10(FDR)") +
    ggplot2::labs(
      x = "ALDEx effect (group contrast)",
      y = "Feature"
    ) +
    ggplot2::theme_bw()
}
