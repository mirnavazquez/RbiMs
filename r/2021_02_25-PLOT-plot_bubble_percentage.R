#' @title Bubble plot of KO/Pathways/Modules and its relative abundance 
#' within each bin.
#' @description Creates a bubble plot of KO/Pathways/Modules relative 
#' abundance within each bin. 
#' It uses the metadata information to color bubbles.
#' @usage plot_bubble_percentage(tibble_ko, x_axis, y_axis, 
#' data_experiment=NULL, color_character=NULL, order_bins=NULL,
#' order_metabolism=NULL, color_pallet=NULL, range_size=NULL)
#' @param tibble_ko a tibble object, created with the mapping_ko 
#' or get_subset_* functions. 
#' @param x_axis a string, a column name of the metabolism table. 
#' It determined the x axis label.
#' @param y_axis a string, a column name of the metabolism table. 
#' It determined the y axis label.
#' @param data_experiment optional. a data frame object 
#' containing metadata information.
#' @param color_character optional. a string column name of the metadata 
#' or metabolism object, used for color.
#' @param order_bins optional. a character vector indicating the bin order.
#' @param order_metabolism optional. a character vector 
#' indicating metabolism order.
#' @param color_pallet optional. a character vector of colors to use.
#' @param range_size optional. a numeric vector indicating 
#' the range size of the dots.
#' @details This function is part of a package used for 
#' the analysis of bins metabolism.
#' @import ggplot2 dplyr rlang pals
#' @examples
#' plot_bubble_percentage(ko_bin_mapp, Bin_name, Module, metadata, Clades)
#' @export
plot_bubble_percentage<-function(tibble_ko,
                                 x_axis,
                                 y_axis,
                                 data_experiment=NULL,
                                 color_character=NULL,
                                 order_bins=NULL,
                                 order_metabolism=NULL,
                                 color_pallet=NULL,
                                 range_size=NULL){
  # Enquoting -------------------------------------------------------------####
  x_axis_enquo <- enquo(x_axis)
  y_axis_enquo <- enquo(y_axis)
  x_axis_label <- as_label(x_axis_enquo)
  y_axis_label <- as_label(y_axis_enquo)
  color_character_enquo <- enquo(color_character)
  # Checking axis ---------------------------------------------------------####
  if( x_axis_label  != "Bin_name") {
    y_axis_enquo<-enquo(x_axis) 
    x_axis_enquo<-enquo(y_axis)
    x_axis_label<-as_label(x_axis_enquo)
    y_axis_label<-as_label(y_axis_enquo)
  } 
  # Checking the color ----------------------------------------------------####
  if(is.null(color_pallet) == T){
    color_pallet<-as.vector(cols25(20))
  }
  # Checking the size -----------------------------------------------------####
  if(is.null(range_size) == T){
    range_size<-c(1,5)
  }
  # Transform from wide to long -------------------------------------------####
  data_to_select<-c("Module", "Module_description", "Pathway", 
                    "Pathway_description", "Genes", 
                    "Gene_description", "Enzyme", "Cycle", 
                    "Pathway_cycle","Detail_cycle", "KO", 
                    "rbims_pathway", "rbims_sub_pathway")
  
  Kegg_long<- tibble_ko %>%
    pivot_longer(cols = -all_of(data_to_select), 
                 values_to = "Abundance",
                 names_to="Bin_name") %>%
    distinct()
  # Count the total number of genes per metabolism-------------------------####
  metabolism_counts_totals<-Kegg_long %>%
    select({{y_axis_enquo}}, .data$KO) %>%
    distinct() %>%
    count(!!y_axis_enquo, sort=T) %>%
    rename(Abundance_metabolism=n)
  # Count the total number of genes per metabolism in genome---------------####
  metabolism_counts_genome<-Kegg_long %>%
    select({{y_axis_enquo}}, .data$Bin_name, .data$KO, 
           .data$Abundance) %>%
    distinct() %>%
    filter(.data$Abundance != "0") %>%
    group_by(.data$Bin_name) %>%
    count(!!y_axis_enquo) %>%
    rename(Abundance_metabolism_genome=n)
  # Calculate the percentage-----------------------------------------------####
  Table_with_percentage<- left_join(
    metabolism_counts_genome, metabolism_counts_totals, 
    by=y_axis_label) %>%
    mutate(Percentage = 
             (.data$Abundance_metabolism_genome * 100) / 
             .data$Abundance_metabolism ) %>%
    select({{y_axis_enquo}}, .data$Bin_name, .data$Percentage) %>%
    distinct() %>%
    drop_na() 
  # Join data experiment --------------------------------------------------####
  if(is.null(data_experiment) == F){
    Table_with_percentage<-Table_with_percentage %>%
      left_join(data_experiment, by="Bin_name")
  }
  # Checking the order Bin ------------------------------------------------####
  if(is.null(order_metabolism) == T){
    order_metabolism<-Table_with_percentage %>%
      ungroup() %>%
      select({{y_axis_enquo}}) %>%
      distinct() %>%
      pull()
  }
  # Checking the order Bin ------------------------------------------------####
  if(is.null(order_bins) == T){
    order_bins<-sort(unique(Table_with_percentage$Bin_name))
  }
  # Checking experiment ---------------------------------------------------####
  if(is.null(data_experiment) == T){
    data_experiment <- NULL
  }
  # Normal / upside -------------------------------------------------------####
  rm(x_axis_enquo, y_axis_enquo, x_axis_label, y_axis_label)
  x_axis_enquo <- enquo(x_axis)
  y_axis_enquo <- enquo(y_axis)
  x_axis_label <- as_label(x_axis_enquo)
  y_axis_label <- as_label(y_axis_enquo)
  # Plot ------------------------------------------------------------------####
  if(x_axis_label == "Bin_name") {
    plot_porcentage<-ggplot(Table_with_percentage,
                            aes(x= factor(!!x_axis_enquo, 
                                          levels = !!order_bins),
                                y= factor(!!y_axis_enquo, 
                                          levels = !!order_metabolism),
                                size= .data$Percentage,
                                color= !!color_character_enquo)) +
      geom_point(alpha=0.5) +
      scale_size(range =range_size) +
      scale_color_manual(values = color_pallet) +
      theme_linedraw() +
      theme(axis.text.x = element_text(size=6, 
                                       angle = 45, 
                                       hjust = 1, 
                                       vjust = 1),
            axis.text.y = element_text(size=5))
  } else if (x_axis_label != "Bin_name" ) {
    plot_porcentage<-ggplot(Table_with_percentage,
                            aes(x= factor(!!x_axis_enquo, 
                                          levels = !!order_metabolism),
                                y= factor(!!y_axis_enquo, 
                                          levels = !!order_bins),
                                size= .data$Percentage,
                                color= !!color_character_enquo)) +
      geom_point(alpha=0.5) +
      scale_size(range = c(1,5)) +
      scale_color_manual(values = color_pallet) +
      theme_linedraw() +
      theme(axis.text.x = element_text(size=6, 
                                       angle = 45, 
                                       hjust = 1, 
                                       vjust = 1),
            axis.text.y = element_text(size=5))
  }
  
  return(plot_porcentage)
}