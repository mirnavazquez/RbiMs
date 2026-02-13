#' @keywords internal
cluster_features <- function(mat_hell){
  # distance matrix
  dist <- vegan::vegdist(mat_hell, method = "bray")
  hcl  <- hclust(dist, method = "ward.D2")
  
  n_feat <- nrow(mat_hell)
  # candidate ks: mínimo 2, máximo 10 o n_feat-1
  ks <- 2:min(10, n_feat - 1)
  
  if (length(ks) == 0) {
    # degenerate case: solo 1 feature, un solo cluster
    k_opt <- 1L
    clusters <- rep(1L, n_feat)
  } else if (length(ks) == 1) {
    k_opt <- ks
    clusters <- cutree(hcl, k_opt)
  } else {
    # compute average silhouette width for each k
    sil_width <- sapply(ks, function(k){
      cl <- cutree(hcl, k)
      mean(cluster::silhouette(cl, dist)[, "sil_width"])
    })
    k_opt <- ks[which.max(sil_width)]
    clusters <- cutree(hcl, k_opt)
  }
  
  list(
    clusters = clusters,
    hclust   = hcl,
    dist     = dist
  )
}
