#' @title Heatmap plot of KO/Pathways/Modules.
#' @description Creates a heatmap of KO/Pathways/Modules of the abundance,
#' presence and absence data or a percentage.
#' @usage heatmap_ko(tibble_ko, y_axis, data_experiment=NULL, calc=NULL,
#' scale_option=NULL, order_y=NULL, order_x=NULL, split_y=FALSE, color_pallet=NULL)
#' @param tibble_ko a tibble object. It could have been
#' created with the mapping_ko or get_subset_* functions. If that is the case
#' you will have to choose which type of calc do you want to show, "binary",
#' "abundance" or "percentage".
#' @param y_axis a string. A column name of the tibble_ko of a 
#' feature to plot (i.e. KO/Pathways/Modules).
#' @param data_experiment a data frame object containing metadata information.
#' @param calc a character indicating with type of calc should 
#' be done to plot the results. Valid values are "Abundance", "Binary", 
#' "Percentage", and "None". If you chose none you are spected to use a
#' tibble table obtained from calc_binary or calc_percentage. 
#' @param scale_option a character indicating if rows or columns should be
#' scale. Valid options "none", "row" or "column".
#' @param order_y is a column name indicating the annotation if rows. 
#' This column name comes from the mapping_ko or get_subset_* object.
#' @param order_x is a column name indicating the annotation if cols. 
#' This column name comes from the metadata object.
#' @param split_y a logical character indicating if you want to split y axis,
#' based on order_y character.
#' @param color_pallet optional. a character vector of colors to use.
#' @details This function is part of a package used for 
#' the calc of bins metabolism.
#' @import pheatmap rlang dplyr tidyr tibble
#' @examples
#' \dontrun{
#' heatmap_ko(tibble_ko=ko_bin_mapp, y_axis=Pathway, 
#' data_experiment = metadata, calc="Percentage", order_y=Module, 
#' order_x=Clades)
#' }
#' @noRd
heatmap_ko<-function(tibble_ko,
                     y_axis,
                     data_experiment=NULL,
                     calc=NULL,
                     scale_option=NULL,
                     order_y=NULL,
                     order_x=NULL,
                     split_y=FALSE,
                     color_pallet=NULL){
  # Enquoting -------------------------------------------------------------####
  y_axis_enquo <- enquo(y_axis)
  order_x_enquo <- enquo(order_x)
  order_y_enquo <- enquo(order_y)
  order_y_label <- as_label(order_y_enquo)
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
  # Check calc --------------------------------------------------------####
  if(calc == "Abundance"){
    tibble_ko_mod<-calc_binary(tibble_ko, !!y_axis_enquo, 
                               binary=FALSE) %>%
      rename(tmp = .data$Abundance)
  } else if (calc == "Binary") {
    tibble_ko_mod<-calc_binary(tibble_ko, !!y_axis_enquo) %>%
      rename(tmp = .data$Presence_absence)
  } else if (calc == "Percentage") {
    tibble_ko_mod<-calc_percentage(tibble_ko, !!y_axis_enquo) %>%
      rename(tmp = .data$Percentage) %>%
      select(-.data$Pathway_number_of_total_elements)
  }  else if (calc == "None") {
    if( "Presence_absence" %in% colnames(tibble_ko)){
      tibble_ko_mod <-rename(tibble_ko, tmp = .data$Presence_absence) %>%
        select({{y_axis_enquo}}, .data$Bin_name, .data$tmp) %>%
        distinct()
      tibble_ko <-rename(tibble_ko, tmp = .data$Presence_absence)
    } else if ( "Abundance" %in% colnames(tibble_ko)){
      tibble_ko_mod <-rename(tibble_ko, tmp = .data$Abundance)%>%
        select({{y_axis_enquo}}, .data$Bin_name, .data$tmp)%>%
        distinct()
      tibble_ko <-rename(tibble_ko, tmp = .data$Abundance)
    } else if ( "Percentage" %in% colnames(tibble_ko)){
      tibble_ko_mod <-rename(tibble_ko, tmp = .data$Percentage)%>%
        select({{y_axis_enquo}}, .data$Bin_name, .data$tmp)%>%
        distinct()
      tibble_ko <-rename(tibble_ko, tmp = .data$Percentage)
    }
  }
  
  # Sorting col names -----------------------------------------------------####
  tibble_ko_colnames<-colnames(tibble_ko)
  paths<-c("Module", "Module_description", "Pathway", 
           "Pathway_description", "Genes", 
           "Gene_description", "Enzyme", "KO", "Cycle", "Pathway_cycle",
           "Detail_cycle", "rbims_pathway", "rbims_sub_pathway")
  if(calc == "None"){
    paths<-c("Module", "Module_description", "Pathway", 
             "Pathway_description", "Genes", 
             "Gene_description", "Enzyme", "KO", "Cycle", "Pathway_cycle",
             "Detail_cycle", "rbims_pathway", "rbims_sub_pathway",
             "tmp", "Bin_name")
  }
  tibble_ko_colnames_2 <- tibble_ko_colnames[tibble_ko_colnames %in%
                                               paths]
  tibble_ko_colnames_3 <- tibble_ko_colnames_2[!tibble_ko_colnames_2 %in%
                                                 y_axis_label]
  # Join ------------------------------------------------------------------####
  table_final<-tibble_ko_mod %>%
    pivot_wider(names_from = .data$Bin_name, 
                values_from = .data$tmp, 
                values_fill = 0) %>%
    distinct() %>%
    left_join(tibble_ko, 
              by= y_axis_label)%>%
    select(-all_of(tibble_ko_colnames_3)) %>%
    select(-contains(".y")) %>%
    distinct() %>%
    rename_all(
      list( ~ stringr::str_replace_all(., ".x", ""))
    ) %>%
    column_to_rownames(y_axis_label)
  # Extracting metabolism -------------------------------------------------####
  if(quo_is_null(order_y_enquo) == F){
    metabolism_order<- tibble_ko %>%
      select({{y_axis_enquo}}, {{order_y_enquo}}) %>%
      drop_na() %>%
      distinct(.data[[y_axis_enquo]], .keep_all = T) %>%
      column_to_rownames(y_axis_label) %>%
      arrange(!!order_y_enquo)
    
  } 
  # Checking the split ----------------------------------------------------####
  if(split_y == FALSE){
    split_y<-NULL
  } else if (split_y == TRUE){
    split_y<-select(metabolism_order, 
                    {{order_y_enquo}})
  }
  # Extracting experiment -------------------------------------------------####
  if(quo_is_null(order_x_enquo) == F && calc != "None"){
    data_to_select<-c("Module", "Module_description", "Pathway", 
                      "Pathway_description", "Genes", 
                      "Gene_description", "Enzyme", "KO", "Cycle", 
                      "Pathway_cycle",
                      "Detail_cycle", "rbims_pathway", "rbims_sub_pathway")
    experiment_order<- tibble_ko %>% 
      pivot_longer(cols = -all_of(data_to_select), 
                   values_to = "Abundance",
                   names_to="Bin_name") %>%
      distinct() %>%
      left_join(data_experiment, by="Bin_name") %>%
      select( .data$Bin_name, {{order_x_enquo}}) %>%
      drop_na() %>%
      distinct(.data$Bin_name, .keep_all = T) %>%
      column_to_rownames("Bin_name") %>%
      arrange(!!order_x_enquo)
  } else if (quo_is_null(order_x_enquo) == F && calc == "None"){
    experiment_order <- tibble_ko %>%
      left_join(data_experiment, by="Bin_name") %>%
      select( .data$Bin_name, {{order_x_enquo}}) %>%
      drop_na() %>%
      distinct(.data$Bin_name, .keep_all = T) %>%
      column_to_rownames("Bin_name") %>%
      arrange(!!order_x_enquo)
  } 
  
  # Order table -----------------------------------------------------------####
  if(quo_is_null(order_x_enquo) == F && quo_is_null(order_y_enquo)== F ){
    sub_samp_ordered <- table_final[rownames(metabolism_order),]
    sub_samp_ordered_2 <- sub_samp_ordered[,rownames(experiment_order)]
  } 
  if(quo_is_null(order_x_enquo) == T && quo_is_null(order_y_enquo)== F ){
    sub_samp_ordered <- table_final[rownames(metabolism_order),]
  }
  if(quo_is_null(order_x_enquo) == F && quo_is_null(order_y_enquo)== T ){
    sub_samp_ordered_2 <- table_final[,rownames(experiment_order)]
  } 
  # Plot ------------------------------------------------------------------####
  if(quo_is_null(order_x_enquo) == F && quo_is_null(order_y_enquo)== F ) {
    plot_heat<-suppressWarnings(
      ComplexHeatmap::pheatmap(sub_samp_ordered_2, 
                               scale = scale_option,
                               annotation_row = metabolism_order,
                               annotation_col = experiment_order,
                               cluster_rows = F,
                               cluster_cols = F,
                               row_split = split_y,
                               column_split = select(experiment_order, 
                                                     {{order_x_enquo}}),
                               main = "Pathway heatmap",
                               fontsize=7,
                               angle_col="0",
                               col = color_pallet)
    )
  }
  # Plot ------------------------------------------------------------------####
  if(quo_is_null(order_x_enquo) == T && quo_is_null(order_y_enquo)== F ) {
    plot_heat<-suppressWarnings(
      ComplexHeatmap::pheatmap(sub_samp_ordered, 
                               scale = scale_option,
                               annotation_row = metabolism_order,
                               cluster_rows = F,
                               cluster_cols = T,
                               row_split = split_y,
                               main = "Pathway heatmap",
                               fontsize=7,
                               angle_col="0",
                               col = color_pallet)
    )
  }
  # Plot ------------------------------------------------------------------####
  if(quo_is_null(order_x_enquo) == F && quo_is_null(order_y_enquo)== T ) {
    plot_heat<-suppressWarnings(
      ComplexHeatmap::pheatmap(sub_samp_ordered_2, 
                               scale = scale_option,
                               annotation_col = experiment_order,
                               cluster_rows = T,
                               cluster_cols = F,
                               column_split = select(experiment_order, 
                                                     {{order_x_enquo}}),
                               main = "Pathway heatmap",
                               fontsize=7,
                               angle_col="45",
                               col = color_pallet)
    )
  }
  # Plot ------------------------------------------------------------------####
  if(quo_is_null(order_x_enquo) == T && quo_is_null(order_y_enquo)== T ) {
    plot_heat<-suppressWarnings(ComplexHeatmap::pheatmap(table_final, 
                                                         scale = NULL,
                                                         cluster_rows = T,
                                                         cluster_cols = T,
                                                         main = "Pathway heatmap",
                                                         fontsize=7,
                                                         angle_col="45",
                                                         col = color_pallet)
    )
    
  }
  
  return(plot_heat)
}


