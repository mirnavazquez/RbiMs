#' @title Bubble plot.
#' @description Creates a bubble plot of KEGG or InterProScan data.
#' @usage plot_bubble(tibble_ko, x_axis, y_axis, analysis=c("KEGG","INTERPRO"),
#' data_experiment=NULL, calc=NULL, color_character=NULL, order_bins=NULL, 
#' order_metabolism=NULL, color_pallet=NULL, range_size=NULL, x_labs=TRUE,
#' y_labs=TRUE, text_x=NULL, text_y=NULL)
#' @param tibble_ko a tibble object, created with the mapping_ko 
#' or get_subset_* functions. 
#' @param x_axis a string, a column name of the metabolism table. 
#' It determined the x axis label.
#' @param y_axis a string, a column name of the metabolism table. 
#' It determined the y axis label.
#' @param analysis a character indicating if your input data are from 
#' KEGG or INTERPRO. 
#' @param calc a character indicating with type of calc should 
#' be done to plot the results. Valid values are "Abundance", "Binary", 
#' "Percentage", and "None". If you chose none you are expected to use a
#' tibble table obtained from calc_binary or calc_percentage. 
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
#' plot_bubble(tibble_ko=ko_bin_mapp, x_axis=Bin_name, y_axis=Module, 
#' analysis="KEGG", data_experiment=metadata, calc="Binary",
#' color_character=Clades)
#' @export
plot_bubble<-function(tibble_ko,
                      x_axis, 
                      y_axis,
                      analysis=c("KEGG","INTERPRO"),
                      data_experiment=NULL,
                      calc=NULL,
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
  # Check analysis --------------------------------------------------------####
  if(analysis == "KEGG"){
    bubble<-bubble_ko(tibble_ko=tibble_ko, x_axis=!!x_axis_enquo, 
                      y_axis=!!y_axis_enquo, calc=calc,
                      data_experiment=data_experiment, 
                      color_character=!!color_character_enquo, 
                      order_bins=order_bins, 
                      order_metabolism=order_metabolism, 
                      color_pallet=color_pallet, range_size=range_size, 
                      x_labs=x_labs, y_labs=y_labs, text_x=text_x, 
                      text_y=text_y)
  }else if (analysis == "INTERPRO") {
    bubble<-bubble_domain(tibble_ko=tibble_ko, x_axis=!!x_axis_enquo, 
                          y_axis=!!y_axis_enquo, data_experiment=data_experiment, 
                          color_character=!!color_character_enquo, 
                          order_bins=order_bins, 
                          order_metabolism=order_metabolism, 
                          color_pallet=color_pallet, range_size=range_size, 
                          x_labs=x_labs,
                          y_labs=y_labs, text_x=text_x, text_y=text_y)
    
  }
  return(bubble)
}