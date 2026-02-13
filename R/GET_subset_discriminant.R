#' Filter feature table by discriminant score
#'
#' @description
#' Wrapper around `get_discriminant_features()` that returns
#' only features with consensus score ≥ threshold.
#' @param tibble_rbims A tibble returned by `consensus_rank()` (or `rbims_consensus()` ...),
#'   containing at least the columns `feature` and `score`.
#' @inheritParams get_discriminant_features
#' @param score_min Minimum consensus score (default = 2)
#' @return Filtered tibble (same format as input) with the full
#' discriminant object attached as attribute "rbims_disc".
#' @export
get_subset_discriminant <- function(tibble_rbims,
                                    analysis,
                                    metadata,
                                    group_col,
                                    score_min    = 1,
                                    min_presence = 2,
                                    feature_col  = NULL) {
  
  # elegir columna de feature si no se especifica
  if (is.null(feature_col)) {
    feature_col <- dplyr::case_when(
      analysis == "KEGG"      ~ "KO",
      analysis == "INTERPRO"  ~ "Pfam",
      analysis == "DBCAN"     ~ "CAZy",
      TRUE                    ~ colnames(tibble_rbims)[1]  # fallback
    )
  }
  
  disc <- get_discriminant_features(
    tibble_profile = tibble_rbims,
    analysis       = analysis,
    feature_col    = feature_col,
    metadata       = metadata,
    group_col      = group_col,
    min_presence   = min_presence
  )
  
  keep <- disc$consensus %>%
    dplyr::filter(score >= score_min) %>%
    dplyr::pull(feature)
  
  out <- tibble_rbims %>%
    dplyr::filter(.data[[feature_col]] %in% keep)
  
  attr(out, "rbims_disc") <- disc
  out
}
