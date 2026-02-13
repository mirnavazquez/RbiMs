#' @keywords internal
differential_features <- function(mat_counts, group){
  stopifnot(ncol(mat_counts) == length(group))
  group <- droplevels(as.factor(group))
  if (nlevels(group) != 2) return(NULL)
  lev <- levels(group)
  pval <- effect <- numeric(nrow(mat_counts))
  for (i in seq_len(nrow(mat_counts))) {
    x <- as.numeric(mat_counts[i, group == lev[1]])
    y <- as.numeric(mat_counts[i, group == lev[2]])
    effect[i] <- mean(log2(y + 1)) - mean(log2(x + 1))
    if (length(unique(c(x, y))) <= 1) {
      pval[i] <- 1
    } else {
      pval[i] <- tryCatch(stats::wilcox.test(x, y)$p.value, error=function(e) 1)
    }
  }
  tibble::tibble(
    feature = rownames(mat_counts),
    effect  = effect,
    p.value = pval,
    we.eBH  = p.adjust(pval, "BH")
  ) |> dplyr::arrange(we.eBH, dplyr::desc(effect))
}
