#' @title Bubble plot of KO/Pathways/Modules abundance within each Bin.
#' @description Creates a bubble plot of KO/Pathways/Modules abundance within each bin. It uses the metadata information to color bubbles.
#' @param tabble_ko a data frame object, created with the mapping_ko or get_subset_* functions.
#' @param other_data a data frame object containing metadata information. 
#' @param x_axis a string, it determined the x axis label.
#' @param y_axis a string, it determined the y axis label.
#' @param size_bubble the Abundance string. It specifies the size of bubbles in the bubble plot. 
#' @param metadata_feature  a string column name of the metadata object, used for color.
#' @details This function is part of a package used for 
#' the analysis of bins metabolism.
#' @import dplyr ggplot2 rlang
#' @examples
#' plot_bubble_abundance(ko_bin_mapp, metadata, Bin_name, Pathway, Abundance, Clades)   
#' @export
plot_bubble_abundance<-function(tabble_ko,
                                other_data,
                                x_axis,
                                y_axis,
                                size_bubble, 
                                metadata_feature){
  ############################ quoting ##############################
  x_axis_enquo <- enquo(x_axis)
  y_axis_enquo <- enquo(y_axis)
  size_bubble_enquo <- enquo(size_bubble)
  metadata_feature_enquo <- enquo(metadata_feature)
  y_axis_label <- as_label(y_axis_enquo)
  
  ######################## Parsing the table ########################
  data_to_select<-c("Module", "Module_description", "Pathway", 
                    "Pathway_description", "Genes", 
                    "Gene_description", "Enzyme", "Cycle", "Pathway_cycle",
                    "Detail_cycle", "KO")
  kegg_longer<-tabble_ko %>%
    pivot_longer(cols = -data_to_select,  
                 names_to = "Bin_name", 
                 values_to = "Abundance") %>%
    distinct() %>%
    mutate_at('Abundance', as.integer) %>%
    select(.data$Bin_name,{{size_bubble_enquo}}, {{y_axis_enquo}} ) %>%
    distinct() %>%
    drop_na() %>%
    left_join(other_data, by="Bin_name") %>%
    mutate(Abundance = case_when(
      .data$Abundance == 0 ~ NA_integer_,
      TRUE ~ as.integer(Abundance)
    )) %>%
    select({{size_bubble_enquo}}, {{metadata_feature_enquo}}, {{y_axis_enquo}}, {{x_axis_enquo}})
  
  ############################## Plot #############################
  Abundance_plot<-ggplot2::ggplot(kegg_longer,
                                  aes(x= !!x_axis_enquo, 
                                      y= !! y_axis_enquo, 
                                      size= !!size_bubble_enquo,
                                      color= !!metadata_feature_enquo)) +
    geom_point(alpha=0.5) +
    scale_size(range = c(1,10))+
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45,
                                     size=8,
                                     hjust = 1,
                                     vjust = 1),
          axis.text.y = element_text(size=7))
  
  return(Abundance_plot)
}