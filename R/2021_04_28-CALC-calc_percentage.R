#' @title Calculate the percentage.
#' @description  Calculate the percentage of KO in certain pathway.
#' @usage calc_percentage(tibble_ko, y_axis, data_experiment=NULL)
#' @param tibble_ko a tibble object from mapping_ko.
#' @param y_axis a character, indicating the pathway to analyze.
#' @param data_experiment optional. a data frame object 
#' containing metadata information.
#' @details This function is part of a package used for the analysis of bins 
#' metabolism.
#' @import  tibble dplyr stringr tidyr  rlang
#' @examples
#' calc_percentage(ko_bin_mapp, Pathway)    
#' @export
calc_percentage<-function(tibble_ko,
                          y_axis,
                          data_experiment=NULL){
  # Enquoting -------------------------------------------------------------####
  y_axis_enquo <- enquo(y_axis)
  y_axis_label <- as_label(y_axis_enquo)
  # Select data -------------------------------------------------------####
  data_to_select<-c("Module", "Module_description", "Pathway", 
                    "Pathway_description", "Genes", 
                    "Gene_description", "Enzyme", "Cycle", "Pathway_cycle",
                    "Detail_cycle", "rbims_pathway", "rbims_sub_pathway", 
                    "KO")
  # Transform from wide to long -------------------------------------------####
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
    rename(Abundance_metabolism=.data$n) %>%
    mutate(new2=str_replace(.data$Abundance_metabolism, "$", "\\)")) %>%
    mutate(new2=str_replace(.data$new2, "^", "\\(")) %>%
    unite("New", {{y_axis_enquo}}, .data$new2, sep = " ", remove=F) %>%
    select(.data$New, {{y_axis_enquo}}, .data$Abundance_metabolism)
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
    select(.data$New, {{y_axis_enquo}}, .data$Bin_name, .data$Percentage) %>%
    distinct() %>%
    drop_na() %>%
    rename( Pathway_number_of_total_elements = .data$New)
  # Join data experiment --------------------------------------------------####
  if(is.null(data_experiment) == F){
    Table_with_percentage<-Table_with_percentage %>%
      left_join(data_experiment, by="Bin_name")
  }
  
  return(Table_with_percentage)
}