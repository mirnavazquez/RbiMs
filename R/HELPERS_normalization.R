#' @keywords internal
normalize_matrices <- function(mat){
  mat_hell <- vegan::decostand(mat, method="hellinger")
  mat_clr  <- compositions::clr(mat + 1)
  mat_prop <- sweep(mat, 2, colSums(mat), "/")
  list(mat_hell = mat_hell, mat_clr = mat_clr, mat_prop = mat_prop)
}
