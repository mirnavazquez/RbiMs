#' Build a KO annotation dictionary from an rbims KEGG table
#'
#' Convenience helper to extract a unique KO → annotation table, including
#' rbims hydrocarbon pathway labels (rbims_pathway, rbims_sub_pathway) when present.
#'
#' @param tibble_rbims A tibble like `ko_hidro_mapp_renamed` containing a KO column
#'   plus optional annotation columns (Genes, Gene_description, etc).
#' @param feature_col Feature ID column name. Default: "KO".
#' @param keep_cols Optional character vector of columns to keep. If NULL, it will
#'   keep a standard set if available.
#'
#' @return A tibble with one row per KO (distinct across kept columns).
#' @export
make_ko_dictionary <- function(tibble_rbims,
                               #' @param keep_cols ...
                               feature_col = "KO",
                               keep_cols = NULL) {
  
  stopifnot(is.data.frame(tibble_rbims))
  stopifnot(feature_col %in% names(tibble_rbims))
  
  # Standard columns we like to keep (only those present will be selected)
  default_cols <- c(
    feature_col,
    "Genes", "Gene_description", "Enzyme",
    "Module", "Module_description",
    "Pathway", "Pathway_description",
    "Cycle", "Pathway_cycle", "Detail_cycle",
    "rbims_pathway", "rbims_sub_pathway"
  )
  
  cols_to_keep <- if (is.null(keep_cols)) default_cols else unique(c(feature_col, keep_cols))
  cols_to_keep <- cols_to_keep[cols_to_keep %in% names(tibble_rbims)]
  
  tibble_rbims %>%
    dplyr::filter(!is.na(.data[[feature_col]])) %>%
    dplyr::select(dplyr::all_of(cols_to_keep)) %>%
    dplyr::distinct()
}






