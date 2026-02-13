#' @keywords internal
indicator_features <- function(mat_hell, clusters){
  # control de permutaciones
  ctrl <- permute::how(nperm = 999)
  
  # multipatt: cada fila de mat_hell = "sitio" (feature), columnas = MAGs
  ind <- indicspecies::multipatt(mat_hell, clusters, control = ctrl)
  
  sign <- ind$sign
  
  # por si cambian algún nombre de columna en versiones futuras
  stat_col <- if ("stat" %in% colnames(sign)) "stat" else colnames(sign)[1]
  p_col    <- if ("p.value" %in% colnames(sign)) "p.value" else colnames(sign)[2]
  
  res <- tibble::tibble(
    feature = rownames(sign),
    indval  = sign[, stat_col],
    p.value = sign[, p_col]
  )
  
  # asignar cluster usando el vector 'clusters'
  res$cluster <- clusters[match(res$feature, rownames(mat_hell))]
  
  res$p.adj <- p.adjust(res$p.value, "BH")
  res
}
