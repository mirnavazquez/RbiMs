#' @title Bubble plot of KO/Pathways/Modules and its relative abundance within each Bin.
#' @description Creates a bubble plot of KO/Pathways/Modules relative abundance within each bin. It uses the metadata information to color bubbles.
#' @param tabble_ko a data frame object, created with the mapping_ko or get_subset_* functions.
#' @param other_data a data frame object containing metadata information. 
#' @param x_axis a string, it determined the x axis label.
#' @param y_axis a string, it determined the y axis label. 
#' @param size_bubble the Percentage string. It specifies the size of bubbles in the bubble plot.  
#' @param metadata_feature a string column name of the metadata object, used for color.
#' @details This function is part of a package used for the analysis of bins metabolism.
#' @import ggplot2 dplyr rlang
#' @examples
#' plot_bubble_percentage(ko_bin_mapp, metadata, Bin_name, Pathway, Percentage, Clades)
#' @export
plot_bubble_percentage<-function(tabble_ko,
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
  #################### Transform from wide to long ####################
  data_to_select<-c("Module", "Module_description", "Pathway", 
                    "Pathway_description", "Genes", 
                    "Gene_description", "Enzyme", "Cycle", "Pathway_cycle",
                    "Detail_cycle", "KO")
  Kegg_long<- tabble_ko %>%
    pivot_longer(cols = -data_to_select, values_to = "Abundance",
                 names_to="Bin_name") %>%
    distinct()
  ######### Count the total number of genes per metabolism ########
  metabolism_counts_totals<-Kegg_long %>%
    select({{y_axis_enquo}}, .data$KO) %>%
    distinct() %>%
    count(!!y_axis_enquo, sort=T) %>%
    rename(Abundance_metabolism=n)
  ## Count the total number of genes per metabolism in each genome ##
  metabolism_counts_genome<-Kegg_long %>%
    select({{y_axis_enquo}}, .data$Bin_name, .data$KO, .data$Abundance) %>%
    distinct() %>%
    filter(.data$Abundance != "0") %>%
    group_by(.data$Bin_name) %>%
    count(!!y_axis_enquo) %>%
    rename(Abundance_metabolism_genome=n)
  ############### Join and calculate the percentage #################
  Table_with_percentage<- left_join(
    metabolism_counts_genome, metabolism_counts_totals, 
    by=y_axis_label) %>%
    mutate(Percentage = (.data$Abundance_metabolism_genome * 100) / .data$Abundance_metabolism ) %>%
    select({{y_axis_enquo}}, .data$Bin_name, .data$Percentage) %>%
    distinct() %>%
    drop_na() %>%
    left_join(other_data, by="Bin_name")
  ################################ Plot #############################
  Table_with_percentage_plot<-ggplot(Table_with_percentage,
                                     aes(x= !!x_axis_enquo, 
                                         y= !! y_axis_enquo, 
                                         size= !!size_bubble_enquo,
                                         color= !!metadata_feature_enquo)) +
    geom_point(alpha=0.5) +
    scale_size(range = c(1,5)) +
    theme_linedraw() +
    theme(axis.text.x = element_text(size=6, 
                                     angle = 45, 
                                     hjust = 1, 
                                     vjust = 1),
          axis.text.y = element_text(size=5))
  
  return(Table_with_percentage_plot)
}
