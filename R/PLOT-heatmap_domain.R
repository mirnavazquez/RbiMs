#' @title Heatmap plot of PFAM/INTERPRO.
#' @description Creates a heatmap of PFAM/INTERPRO of the abundance,
#' presence and absence data or a percentage.
#' @usage heatmap_domain(tibble_ko, y_axis, scale_option=NULL, 
#' color_pallet=NULL, distance=FALSE)
#' @param tibble_ko a tibble object. It could have been
#' created with the read_interpro or get_subset_* functions.
#' @param y_axis a string. A column name of the tibble_ko of a 
#' feature to plot (i.e. PFAM/INTERPRO).
#' @param scale_option a character indicating if rows or columns should be
#' scale. Valid options "none", "row" or "column".
#' @param color_pallet optional. a character vector of colors to use.
#' @param distance optional. If TRUE it will calculate a distance matrix and 
#' show how similar are the different genomes.
#' @details This function is part of a package used for 
#' the analysis of bins metabolism.
#' @import pheatmap rlang dplyr tidyr tibble
#' @examples
#' \dontrun{
#' heatmap_domain(input_data_profile, y_axis=PFAM, 
#' scale_option="none", distance=T)
#' }
#' @noRd
heatmap_domain<-function(tibble_ko, 
                         y_axis,
                         scale_option=NULL, 
                         color_pallet=NULL,
                         distance=FALSE){
  # Enquoting -------------------------------------------------------------####
  y_axis_enquo <- enquo(y_axis)
  y_axis_label <- as_label(y_axis_enquo)
  # Checking the scale ----------------------------------------------------####
  if(is.null(scale_option) == T){
    scale_option<-"none"
  } else if (is.null(scale_option) == "row"){
    scale_option<-"row"
  }else if (is.null(scale_option) == "column"){
    scale_option<-"column"
  }
  # Checking the color ----------------------------------------------------####
  if(is.null(color_pallet) == T){
    color_pallet<-viridis(n=100)
  }
  # Preparing dataframe ---------------------------------------------------####
  heatmap_domain_table<-tibble_ko %>%
    select(-.data$domain_name) %>%
    column_to_rownames({{y_axis_label}})
  # Checking distance -----------------------------------------------------####
  if(isTRUE(distance) == T) {
    heatmap_domain_table_2<-tibble_ko %>%
      select(-.data$domain_name) %>%
      pivot_longer(cols = -{{y_axis_enquo}}, names_to="Bin_name", 
                   values_to = "count") %>%
      pivot_wider(names_from = {{y_axis_enquo}}, values_from = count, 
                  values_fill = 0) %>%
      column_to_rownames("Bin_name")
    distance_domains<-stats::dist(heatmap_domain_table_2)
    heatmap_domain_table<-as.matrix(distance_domains, method="euclidean")
  }
  # Plot ------------------------------------------------------------------####
  plot_heat<-pheatmap::pheatmap(heatmap_domain_table, 
                                angle_col="45", 
                                scale = scale_option,
                                main = "Protein domain heatmap",
                                col = color_pallet,
                                cluster_rows = T,
                                cluster_cols = T)
  
  return(plot_heat)
}