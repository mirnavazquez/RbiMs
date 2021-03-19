#' @title Bubble plot of KO/Pathways/Modules and its relative abundance within each Bin.
#' @description Creates a bubble plot of KO/Pathways/Modules relative abundance within each bin. It uses the metadata information to color bubbles.
#' @param tabble_ko a data frame object, created with the mapping_ko or get_subset_* functions.
#' @param other_data a data frame object containing metadata information. 
#' @param x_axis a string, it determined the x axis label.
#' @param y_axis a string, it determined the y axis label. 
#' @param metadata_feature a string column name of the metadata object, used for color.
#' @param order_bins a character vector indicating the bin order.
#' @param order_metabolism a character vector indicaton metabolism order.
#' @param color_bin a character vector of colors to use.
#' @details This function is part of a package used for the analysis of bins metabolism.
#' @import ggplot2 dplyr rlang pals
#' @examples
#' plot_bubble_percentage(ko_bin_mapp, metadata, Bin_name, Module, Clades)
#' @export
plot_bubble_percentage<-function(tabble_ko,
                                 other_data,
                                 x_axis,
                                 y_axis,
                                 metadata_feature=NULL,
                                 order_bins=NULL,
                                 order_metabolism=NULL,
                                 color_bin=NULL){
  #################################################################
  ######################## Enquoting ##############################
  #################################################################
  x_axis_enquo <- enquo(x_axis)
  y_axis_enquo <- enquo(y_axis)
  
  x_axis_label <- as_label(x_axis_enquo)
  y_axis_label <- as_label(y_axis_enquo)
  
  if( x_axis_label  != "Bin_name") {
    
    y_axis_enquo<-enquo(x_axis) 
    x_axis_enquo<-enquo(y_axis)
    
    x_axis_label<-as_label(x_axis_enquo)
    y_axis_label<-as_label(y_axis_enquo)
    
  } 
  
  
  metadata_feature_enquo <- enquo(metadata_feature)
  
  if(is.null(color_bin) == T){
    color_bin<-as.vector(cols25(20))
  }
  
  #################################################################
  ####### Transform from wide to long #############################
  #################################################################
  data_to_select<-c("Module", "Module_description", "Pathway", 
                    "Pathway_description", "Genes", 
                    "Gene_description", "Enzyme", "Cycle", 
                    "Pathway_cycle","Detail_cycle", "KO", 
                    "rbims_pathway", "rbims_sub_pathway")
  
  Kegg_long<- tabble_ko %>%
    pivot_longer(cols = -all_of(data_to_select), 
                 values_to = "Abundance",
                 names_to="Bin_name") %>%
    distinct()
  #################################################################
  ####### Count the total number of genes per metabolism ##########
  #################################################################
  metabolism_counts_totals<-Kegg_long %>%
    select({{y_axis_enquo}}, .data$KO) %>%
    distinct() %>%
    count(!!y_axis_enquo, sort=T) %>%
    rename(Abundance_metabolism=n)
  #################################################################
  #Count the total number of genes per metabolism in each genome ##
  #################################################################
  metabolism_counts_genome<-Kegg_long %>%
    select({{y_axis_enquo}}, .data$Bin_name, .data$KO, 
           .data$Abundance) %>%
    distinct() %>%
    filter(.data$Abundance != "0") %>%
    group_by(.data$Bin_name) %>%
    count(!!y_axis_enquo) %>%
    rename(Abundance_metabolism_genome=n)
  #################################################################
  ############### Join and calculate the percentage ###############
  #################################################################
  Table_with_percentage<- left_join(
    metabolism_counts_genome, metabolism_counts_totals, 
    by=y_axis_label) %>%
    mutate(Percentage = 
             (.data$Abundance_metabolism_genome * 100) / 
             .data$Abundance_metabolism ) %>%
    select({{y_axis_enquo}}, .data$Bin_name, .data$Percentage) %>%
    distinct() %>%
    drop_na() %>%
    left_join(other_data, by="Bin_name")
  #################################################################
  ####################### Pot Normal / upsite #####################
  #################################################################
  rm(x_axis_enquo, y_axis_enquo, x_axis_label, y_axis_label)
  x_axis_enquo <- enquo(x_axis)
  y_axis_enquo <- enquo(y_axis)
  
  x_axis_label <- as_label(x_axis_enquo)
  y_axis_label <- as_label(y_axis_enquo)
  
  
  if(x_axis_label == "Bin_name") {
    if(is.null(order_bins) == T && is.null(order_metabolism) == T ){
      Table_with_percentage_plot<-ggplot(Table_with_percentage,
                                         aes(x= !!x_axis_enquo, 
                                             y= !!y_axis_enquo, 
                                             size= .data$Percentage,
                                             color= !!metadata_feature_enquo)) +
        geom_point(alpha=0.5) +
        scale_size(range = c(1,5)) +
        scale_color_manual(values = color_bin) +
        theme_linedraw() +
        theme(axis.text.x = element_text(size=6, 
                                         angle = 45, 
                                         hjust = 1, 
                                         vjust = 1),
              axis.text.y = element_text(size=5))
    } else if (is.null(order_bins) == F && is.null(order_metabolism) == T) {
      Table_with_percentage_plot<-ggplot(Table_with_percentage,
                                         aes(x= factor(!!x_axis_enquo, 
                                                       levels = order_bins),
                                             y= !!y_axis_enquo, 
                                             size= .data$Percentage,
                                             color= !!metadata_feature_enquo)) +
        geom_point(alpha=0.5) +
        scale_size(range = c(1,5)) +
        scale_color_manual(values = color_bin) +
        theme_linedraw() +
        theme(axis.text.x = element_text(size=6, 
                                         angle = 45, 
                                         hjust = 1, 
                                         vjust = 1),
              axis.text.y = element_text(size=5))   
    } else if(is.null(order_bins) == F && is.null(order_metabolism) == F) {
      Table_with_percentage_plot<-ggplot(Table_with_percentage,
                                         aes(x= factor(!!x_axis_enquo, levels = !!order_bins),
                                             y= factor(!!y_axis_enquo, levels = !!order_metabolism),
                                             size= .data$Percentage,
                                             color= !!metadata_feature_enquo)) +
        geom_point(alpha=0.5) +
        scale_size(range = c(1,5)) +
        scale_color_manual(values = color_bin) +
        theme_linedraw() +
        theme(axis.text.x = element_text(size=6, 
                                         angle = 45, 
                                         hjust = 1, 
                                         vjust = 1),
              axis.text.y = element_text(size=5))
    } else if(is.null(order_bins) == T && is.null(order_metabolism) == F){
      Table_with_percentage_plot<-ggplot(Table_with_percentage,
                                         aes(x= !!x_axis_enquo,
                                             y= factor(!!y_axis_enquo, levels = !!order_metabolism),
                                             size= .data$Percentage,
                                             color= !!metadata_feature_enquo)) +
        geom_point(alpha=0.5) +
        scale_size(range = c(1,5)) +
        scale_color_manual(values = color_bin) +
        theme_linedraw() +
        theme(axis.text.x = element_text(size=6, 
                                         angle = 45, 
                                         hjust = 1, 
                                         vjust = 1),
              axis.text.y = element_text(size=5))
    }
    
  } else if (x_axis_label != "Bin_name" ) {
    if(is.null(order_bins) == T && is.null(order_metabolism) == T ){
      Table_with_percentage_plot<-ggplot(Table_with_percentage,
                                         aes(x= !!x_axis_enquo, 
                                             y= !!y_axis_enquo, 
                                             size= .data$Percentage,
                                             color= !!metadata_feature_enquo)) +
        geom_point(alpha=0.5) +
        scale_size(range = c(1,5)) +
        scale_color_manual(values = color_bin) +
        theme_linedraw() +
        theme(axis.text.x = element_text(size=6, 
                                         angle = 45, 
                                         hjust = 1, 
                                         vjust = 1),
              axis.text.y = element_text(size=5))
    } else if (is.null(order_bins) == F && is.null(order_metabolism) == T) {
      Table_with_percentage_plot<-ggplot(Table_with_percentage,
                                         aes(x= factor(!!x_axis_enquo, levels = !!order_bins),
                                             y= !!y_axis_enquo, 
                                             size= .data$Percentage,
                                             color= !!metadata_feature_enquo)) +
        geom_point(alpha=0.5) +
        scale_size(range = c(1,5)) +
        scale_color_manual(values = color_bin) +
        theme_linedraw() +
        theme(axis.text.x = element_text(size=6, 
                                         angle = 45, 
                                         hjust = 1, 
                                         vjust = 1),
              axis.text.y = element_text(size=5))    
    } else if(is.null(order_bins) == F && is.null(order_metabolism) == F) {
      Table_with_percentage_plot<-ggplot(Table_with_percentage,
                                         aes(x= factor(!!x_axis_enquo, levels = !!order_bins),
                                             y= factor(!!y_axis_enquo, levels = !!order_metabolism),
                                             size= .data$Percentage,
                                             color= !!metadata_feature_enquo)) +
        geom_point(alpha=0.5) +
        scale_size(range = c(1,5)) +
        scale_color_manual(values = color_bin) +
        theme_linedraw() +
        theme(axis.text.x = element_text(size=6, 
                                         angle = 45, 
                                         hjust = 1, 
                                         vjust = 1),
              axis.text.y = element_text(size=5))
    } else if(is.null(order_bins) == T && is.null(order_metabolism) == F){
      Table_with_percentage_plot<-ggplot(Table_with_percentage,
                                         aes(x= !!x_axis_enquo,
                                             y= factor(!!y_axis_enquo, levels = !!order_metabolism),
                                             size= .data$Percentage,
                                             color= !!metadata_feature_enquo)) +
        geom_point(alpha=0.5) +
        scale_size(range = c(1,5)) +
        scale_color_manual(values = color_bin) +
        theme_linedraw() +
        theme(axis.text.x = element_text(size=6, 
                                         angle = 45, 
                                         hjust = 1, 
                                         vjust = 1),
              axis.text.y = element_text(size=5))
    }
  }
  
  return(Table_with_percentage_plot)
}