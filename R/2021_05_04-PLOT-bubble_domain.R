#' @title Bubble plot of PFAM/INTERPRO and its relative abundance 
#' within each bin.
#' @description Creates a bubble plot of PFAM/INTERPRO relative 
#' abundance within each bin. 
#' It uses the metadata information to color bubbles.
#' @usage bubble_domain(tibble_ko, x_axis, y_axis, data_experiment=NULL,
#' color_character=NULL, order_bins=NULL, order_metabolism=NULL, 
#' color_pallet=NULL, range_size=NULL, x_labs=TRUE, y_labs=TRUE,
#' text_x=NULL, text_y=NULL)
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
#' @param x_labs optional. If FALSE it will set the x lab to NULL. 
#' @param y_labs optional. If FALSE it will set the y lab to NULL. 
#' @param text_x optional. A numeric vector indicating the size
#'  of the x text letters.
#' @param text_y optional. A numeric vector indicating the size
#'  of the y text letters.
#' @details This function is part of a package used for 
#' the analysis of bins metabolism.
#' @import ggplot2 dplyr rlang pals
#' @examples
#' \dontrun{
#' bubble_domain(tibble_ko=input_data, x_axis=Bin_name, y_axis=PFAM, 
#' data_experiment=metadata, color_character=Genus)
#' }
#' @noRd
bubble_domain<-function(tibble_ko,
                        x_axis, 
                        y_axis,
                        data_experiment=NULL,
                        color_character=NULL,
                        order_bins=NULL,
                        order_metabolism=NULL,
                        color_pallet=NULL, 
                        range_size=NULL,
                        x_labs=TRUE,
                        y_labs=TRUE,
                        text_x=NULL,
                        text_y=NULL){
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
    range_size<-c(1,10)
  }
  # Checking the xlabs ----------------------------------------------------####
  if(isTRUE(x_labs) == T){
    x_labs<-x_axis_enquo
  } else if (isTRUE(x_labs) == F) {
    x_labs<-NULL
  }
  # Checking the ylabs ----------------------------------------------------####
  if(isTRUE(y_labs) == T){
    y_labs<-y_axis_enquo
  } else if (isTRUE(y_labs) == F) {
    y_labs<-NULL
  }
  # Checking the text_x ----------------------------------------------------####
  if(is.null(text_x) == T){
    text_x<-7
  }
  # Checking the text_x ----------------------------------------------------####
  if(is.null(text_y) == T){
    text_y<-7
  }
  # Table -----------------------------------------------------------------####
  bubble<-tibble_ko %>%
    select(-.data$domain_name) %>%
    pivot_longer(cols = -{{y_axis_enquo}}, names_to="Bin_name", 
                 values_to = "Abundance") %>%
    mutate_at('Abundance', as.integer) %>%
    mutate(Abundance = case_when(
      .data$Abundance == 0 ~ NA_integer_,
      TRUE ~ as.integer(Abundance)
    )) 
  # Join data experiment --------------------------------------------------####
  if(is.null(data_experiment) == F){
    bubble<-bubble %>%
      left_join(data_experiment, by="Bin_name")
  }
  # Checking the order ---------------------------------------------------####
  if(is.null(order_metabolism) == T){
    order_metabolism<-bubble %>%
      ungroup() %>%
      select({{y_axis_enquo}}) %>%
      distinct() %>%
      pull()
  }
  # Checking the order ---------------------------------------------------####
  if(is.null(order_bins) == T){
    order_bins<-sort(unique(bubble$Bin_name))
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
    plot_bubble<-ggplot(bubble,
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
      theme(axis.text.x = element_text(size=text_x, 
                                       angle = 45, 
                                       hjust = 1, 
                                       vjust = 1),
            axis.text.y = element_text(size=text_y))+
      xlab(x_labs) + 
      ylab(y_labs)
  } else if (x_axis_label != "Bin_name" ) {
    plot_bubble<-ggplot(bubble,
                        aes(x= factor(!!x_axis_enquo, 
                                      levels = !!order_metabolism),
                            y= factor(!!y_axis_enquo, 
                                      levels = !!order_bins),
                            size= .data$Abundance,
                            color= !!color_character_enquo)) +
      geom_point(alpha=0.5) +
      scale_size(range = range_size) +
      scale_color_manual(values = color_pallet) +
      theme_linedraw() +
      theme(axis.text.x = element_text(size=text_x, 
                                       angle = 45, 
                                       hjust = 1, 
                                       vjust = 1),
            axis.text.y = element_text(size=text_y)) +
      xlab(x_labs) + 
      ylab(y_labs) 
  }
  
  return(plot_bubble)
}