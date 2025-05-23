#' @title Heatmap plot of MEROPS.
#' @description Creates a heatmap of the abundance of MEROPS families,
#' presence and absence data or a percentage.
#' @usage heatmap_merops(tibble_ko, y_axis, scale_option=NULL, 
#' color_pallet=NULL, distance=FALSE)
#' @param tibble_ko a tibble object. It could have been
#' created with the read_interpro or get_subset_* functions.
#' @param y_axis a string. A column name of the tibble_ko of a 
#' feature to plot (i.e. MEROPS_family/domain_name).
#' @param scale_option a character indicating if rows or columns should be
#' scaled. Valid options "none", "row" or "column".
#' @param color_pallet optional. A character vector of colors to use.
#' @param distance optional. If TRUE, it will calculate a distance matrix and 
#' show how similar are the different genomes.
#' @details This function is part of a package used for 
#' the analysis of MEROPS profiles.
#' @import ComplexHeatmap rlang dplyr tidyr tibble
#' @examples
#' \dontrun{
#' heatmap_merops(input_data_profile, y_axis="domain_name", 
#' scale_option="none", distance=T)
#' }
#' @noRd
heatmap_merops <- function(tibble_ko, 
                           y_axis,
                           scale_option = NULL, 
                           color_pallet = NULL,
                           distance = FALSE) {
  
  # Enquoting -------------------------------------------------------------####
  y_axis_enquo <- enquo(y_axis)
  y_axis_label <- as_label(y_axis_enquo)
  
  # Checking the scale ----------------------------------------------------####
  if (is.null(scale_option)) {
    scale_option <- "none"
  } else if (scale_option == "row") {
    scale_option <- "row"
  } else if (scale_option == "column") {
    scale_option <- "column"
  }
  
  # Checking the color ----------------------------------------------------####
  if (is.null(color_pallet)) {
    color_pallet <- viridis::viridis(n = 100)
  }
  
  # Preparing dataframe ---------------------------------------------------####
  heatmap_merops_table <- tibble_ko %>%
    select(-domain_name) %>%
    column_to_rownames({{y_axis_label}})
  
  # Checking distance -----------------------------------------------------####
  if (isTRUE(distance)) {
    heatmap_merops_table_2 <- tibble_ko %>%
      select(-domain_name) %>%
      pivot_longer(cols = -{{y_axis_enquo}}, names_to = "Bin_name", 
                   values_to = "count") %>%
      pivot_wider(names_from = {{y_axis_enquo}}, values_from = count, 
                  values_fill = 0) %>%
      column_to_rownames("Bin_name")
    
    distance_domains <- stats::dist(heatmap_merops_table_2)
    heatmap_merops_table <- as.matrix(distance_domains, method = "euclidean")
  }
  
  # Plot ------------------------------------------------------------------####
  plot_heat <- ComplexHeatmap::pheatmap(heatmap_merops_table, 
                                        angle_col = 45, 
                                        name = "Abundance",
                                        scale = scale_option,
                                        main = "MEROPS Families Heatmap",
                                        col = color_pallet,
                                        cluster_rows = TRUE,
                                        cluster_cols = TRUE)
  
  return(plot_heat)
}
