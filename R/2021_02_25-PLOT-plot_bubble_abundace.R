#' @title Bubble plot of KO/Pathways/Modules abundance within each Bin.
#' @description Creates a bubble plot of KO/Pathways/Modules abundance 
#' within each bin. It uses the metadata information to color bubbles.
#' @param tibble_ko a data frame object, created with 
#' the mapping_ko or get_subset_* functions.
#' @usage plot_bubble_abundance(tibble_ko, x_axis, y_axis, 
#' data_experiment=NULL, color_character=NULL, 
#' order_bins=NULL, order_metabolism=NULL, 
#' color_pallet=NULL, range_size=NULL)
#' @param x_axis a string, it determined the x axis label.
#' @param y_axis a string, it determined the y axis label.
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
#' @import dplyr ggplot2 rlang
#' @examples
#' plot_bubble_abundance(ko_bin_mapp, Bin_name, Pathway)    
#' @export
plot_bubble_abundance<-function(tibble_ko,
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
  # Calculate abundance ---------------------------------------------------####
  data_to_select<-c("Module", "Module_description", "Pathway", 
                    "Pathway_description", "Genes", 
                    "Gene_description", "Enzyme", "Cycle", "Pathway_cycle",
                    "Detail_cycle", "KO", "rbims_pathway", "rbims_sub_pathway")
  kegg_longer<-tibble_ko %>%
    pivot_longer(cols = -all_of(data_to_select),  
                 names_to = "Bin_name", 
                 values_to = "Abundance") %>%
    distinct() %>%
    mutate_at('Abundance', as.integer) %>%
    select(.data$Bin_name, .data$Abundance, {{y_axis_enquo}} ) %>%
    distinct() %>%
    drop_na()
  # Evaluate experiment ---------------------------------------------------####
  if(is.null(data_experiment) == F){
    kegg_longer<-left_join(kegg_longer, data_experiment, by="Bin_name")
  }
  # Calculate abundance ---------------------------------------------------####
  kegg_longer<-kegg_longer %>%
    mutate(Abundance = case_when(
      .data$Abundance == 0 ~ NA_integer_,
      TRUE ~ as.integer(Abundance)
    )) %>%
    select(.data$Abundance, {{color_character_enquo}}, 
           {{y_axis_enquo}}, {{x_axis_enquo}})
  # Checking the order ---------------------------------------------------####
  if(is.null(order_metabolism) == T){
    order_metabolism<-kegg_longer %>%
      ungroup() %>%
      select({{y_axis_enquo}}) %>%
      distinct() %>%
      pull()
  }
  # Checking the order ---------------------------------------------------####
  if(is.null(order_bins) == T){
    order_bins<-sort(unique(kegg_longer$Bin_name))
  }
  # Plot ------------------------------------------------------------------####
    if(x_axis_label == "Bin_name") {
      Abundance_plot<-ggplot(kegg_longer,
                             aes(x= factor(!!x_axis_enquo,
                                           levels = !!order_bins),
                                 y= factor(!!y_axis_enquo,
                                           levels = !!order_metabolism),
                                 size= .data$Abundance,
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
      Abundance_plot<-ggplot(kegg_longer,
                             aes(x= factor(!!x_axis_enquo,
                                           levels = !!order_bins),
                                 y= factor(!!y_axis_enquo,
                                           levels = !!order_metabolism),
                                 size= .data$Abundance,
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
    }
  
  return(Abundance_plot)
}


