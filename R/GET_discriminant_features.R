#' Identify discriminant features across MAGs
#'
#' @description
#' Integrates clustering, indicator analysis, differential testing,
#' and Random Forest importance into a unified workflow.
#'
#' @param tibble_profile Wide profile table (features × MAGs)
#' @param analysis Character: "KEGG", "Pfam", "INTERPRO", "dbCAN", or "MEROPS"
#' @param feature_col Column name containing feature IDs
#' @param metadata MAG-level metadata (rownames = MAG identifiers)
#' @param group_col Metadata column to discriminate (e.g. "Depth", "Class", "Phylum")
#' @param norm Normalization method ("hellinger" or "clr")
#' @param min_presence Minimum number of MAGs in which a feature must appear
#' @return List with: matrices (counts, hellinger, clr), clustering,
#' indicator results, differential tests, Random Forest importance, and consensus ranking.
#' @export
get_discriminant_features <- function(
    tibble_profile,
    analysis,
    feature_col,
    metadata,
    group_col = NULL,
    norm = "hellinger",
    min_presence = 2
){
  # --- Feature matrix ---
  mat <- rbims_feature_matrix(tibble_profile, feature_col)
  mat_filtered <- mat[rowSums(mat > 0) >= min_presence, ]
  
  # --- Normalization ---
  norm_mats <- normalize_matrices(mat_filtered)
  
  # --- Clustering ---
  clust <- cluster_features(norm_mats$mat_hell)
  
  # --- Indicator analysis ---
  ind <- indicator_features(norm_mats$mat_hell, clust$clusters)
  
  # --- Differential analysis ---
  diff <- NULL
  if (!is.null(group_col) && group_col %in% colnames(metadata)) {
    diff <- differential_features(mat_filtered, metadata[[group_col]])
  }
  
  # --- Random Forest importance ---
  rf <- NULL
  if (!is.null(group_col) && group_col %in% colnames(metadata)) {
    rf <- rf_importance_features(norm_mats$mat_prop, metadata[[group_col]])
  }
  
  # --- Consensus ---
  consensus <- consensus_rank(ind, diff, rf)
  
  list(
    analysis       = analysis,
    feature_col    = feature_col,
    mat_counts     = mat_filtered,
    mat_hell       = norm_mats$mat_hell,
    mat_clr        = norm_mats$mat_clr,
    clusters       = clust$clusters,
    hclust         = clust$hclust,
    dist           = clust$dist,
    indicators     = ind,
    differential   = diff,
    rf_importance  = rf,
    consensus      = consensus
  )
}
