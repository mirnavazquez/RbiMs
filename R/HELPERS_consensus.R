#' @keywords internal
consensus_rank <- function(ind_df, diff_df = NULL, rf_df = NULL){
  # 1) lista de features presentes en cualquiera de las tablas
  all_features <- unique(c(
    ind_df$feature,
    if (!is.null(diff_df)) diff_df$feature,
    if (!is.null(rf_df))   rf_df$feature
  ))
  
  out <- tibble::tibble(feature = all_features)
  
  # 2) indicadores (siempre existen como data.frame, aunque sin señal)
  out <- out |>
    dplyr::left_join(
      ind_df |> dplyr::select(feature, ind_padj = p.adj),
      by = "feature"
    )
  
  # 3) análisis diferencial (opcional)
  if (!is.null(diff_df)) {
    out <- out |>
      dplyr::left_join(
        diff_df |> dplyr::select(feature, effect, aldex_FDR = we.eBH),
        by = "feature"
      )
  } else {
    out$effect    <- NA_real_
    out$aldex_FDR <- NA_real_
  }
  
  # 4) Random Forest (opcional)
  if (!is.null(rf_df)) {
    out <- out |>
      dplyr::left_join(
        rf_df |> dplyr::select(feature, rf_importance),
        by = "feature"
      )
  } else {
    out$rf_importance <- NA_real_
  }
  
  # 5) umbral de importancia RF (solo si hay valores no-NA)
  if (all(is.na(out$rf_importance))) {
    rf_thr <- Inf
  } else {
    rf_thr <- stats::quantile(out$rf_importance, 0.75, na.rm = TRUE)
  }
  
  # 6) flags (tratando NA como FALSE) + score
  out <- out |>
    dplyr::mutate(
      flag_ind  = !is.na(ind_padj)   & ind_padj  < 0.05,
      flag_diff = !is.na(aldex_FDR)  & aldex_FDR < 0.05,
      flag_rf   = !is.na(rf_importance) & rf_importance > rf_thr,
      score     = as.integer(flag_ind) +
        as.integer(flag_diff) +
        as.integer(flag_rf)
    ) |>
    dplyr::arrange(
      dplyr::desc(score),
      dplyr::desc(abs(dplyr::coalesce(effect, 0))),
      dplyr::desc(dplyr::coalesce(rf_importance, 0))
    )
  
  out
}

