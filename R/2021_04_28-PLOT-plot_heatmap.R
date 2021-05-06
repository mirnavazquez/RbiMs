#' @title Heatmap
#' @description Creates a heatmap of KEGG and InterProScan data.
#' @usage plot_heatmap(tibble_ko, y_axis, analysis=c("KEGG","INTERPRO"), 
#' data_experiment=NULL, calc=NULL, scale_option=NULL, order_y=NULL, order_x=NULL,
#' split_y=FALSE, color_pallet=NULL, distance=FALSE)
#' @param tibble_ko a data frame object. It could have been
#' created with the mapping_ko or get_subset_* functions. If that is the case
#' you will have to choose which type of analysis do you want to show, "binary",
#' "abundance" or "percentage". You can also include a table with the counts that 
#' you want to use.
#' @param y_axis a string. A column name of the tibble_ko of a 
#' feature to plot (i.e. KO/Pathways/Modules/PFAM/INTERPRO).
#' @param analysis a character indicating the input file. Valid arguments "KEGG" or
#' "INTERPRO".
#' @param data_experiment optional. a data frame object containing metadata information.
#' @param calc optional. This option is only valid for "KEGG" analysis. 
#' Is a character indicating with type of calc should 
#' be done to plot the results. Valid values are "Abundance", "Binary", 
#' "Percentage", and "None". If you chose none you are spected to use a
#' tibble table obtained from calc_binary or calc_percentage.
#' @param scale_option a character indicating if rows or columns should be
#' scale. Valid options "none", "row" or "column".
#' @param order_y optional. This option is only valid for "KEGG" analysis.
#' Is a column name indicating the annotation if rows. 
#' This column name comes from the mapping_ko or get_subset_* object.
#' @param order_x optional. This option is only valid for "KEGG" analysis.
#' Is a column name indicating the annotation if cols. 
#' This column name comes from the metadata object.
#' @param split_y optional. This option is only valid for "KEGG" analysis. 
#' A logical character indicating if you want to split y axis,
#' based on order_y character.
#' @param color_pallet optional. a character vector of colors to use.
#' @param distance optional. If TRUE it will calculate a distance matrix and 
#' show how similar are the different genomes.
#' @details This function is part of a package used for 
#' the analysis of bins metabolism.
#' @import pheatmap rlang dplyr tidyr tibble
#' @examples
#' \dontrun{
#' plot_heatmap(tibble_ko=ko_bin_mapp, y_axis=Pathway,
#' analysis="INTERPRO", data_experiment = metadata, 
#' calc="Binary", order_y=Module, order_x=Clades)
#' }
#' @export
plot_heatmap<-function(tibble_ko,
                       y_axis,
                       analysis=c("KEGG","INTERPRO"),
                       data_experiment=NULL,
                       calc=NULL,
                       scale_option=NULL,
                       order_y=NULL,
                       order_x=NULL,
                       split_y=FALSE,
                       color_pallet=NULL,
                       distance=FALSE){
  # Enquoting -------------------------------------------------------------####
  y_axis_enquo <- enquo(y_axis)
  order_x_enquo <- enquo(order_x)
  order_y_enquo <- enquo(order_y)
  order_y_label <- as_label(order_y_enquo)
  y_axis_label <- as_label(y_axis_enquo)
  # Check analysis --------------------------------------------------------####
  if(analysis == "KEGG"){
    plot_heat<-heatmap_ko(tibble_ko, !!y_axis_enquo, 
                          data_experiment=data_experiment, 
                          calc=calc,scale_option=scale_option, 
                          order_y=!!order_y_enquo, 
                          order_x=!!order_x_enquo, split_y=split_y, 
                          color_pallet=color_pallet)
  }else if (analysis == "INTERPRO") {
    plot_heat<-heatmap_domain(tibble_ko, !!y_axis_enquo, 
                              scale_option=scale_option, 
                              color_pallet=color_pallet, distance=distance)
    
  }
  return(plot_heat)
}