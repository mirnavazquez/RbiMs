#' @keywords internal
rbims_feature_matrix <- function(tibble_profile, feature_col){
  # solo quedarnos con la columna de feature + columnas numéricas
  mat_df <- tibble_profile %>%
    dplyr::select(dplyr::all_of(feature_col), where(is.numeric)) %>%
    dplyr::group_by(.data[[feature_col]]) %>%
    dplyr::summarise(dplyr::across(where(is.numeric), sum, na.rm = TRUE),
                     .groups = "drop")
  
  mat <- mat_df %>%
    tibble::column_to_rownames(feature_col) %>%
    as.matrix()
  
  storage.mode(mat) <- "numeric"
  mat
}

