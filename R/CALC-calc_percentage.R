#' @title Calculate the percentage.
#' @description Calculate the percentage of a certain pathway in the data. 
#' @usage calc_percentage (tibble_ko, y_axis, data_experiment = NULL)
#' @param tibble_ko a tibble object from mapping_ko.
#' @param y_axis a character, indicating the pathway to analyze.
#' @param data_experiment optional. a data frame object containing metadata 
#' information.
#' @details Calculate the percentage of a certain pathway in the data.
#' @import  tibble dplyr stringr tidyr  rlang
#' @examples 
#' # calc_percentage (ko_bin_mapp, Pathway)   
#' @return A data frame with the calculated percentage of genes in a certain 
#' pathway of the data.
#' @export

#Sets the syntax of the function-------------------------------------------####
calc_percentage <- function(tibble_ko,
                            y_axis,
                            data_experiment = NULL){
# Enquoting ---------------------------------------------------------------####
  y_axis_enquo <- enquo(y_axis)
  y_axis_label <- as_label(y_axis_enquo)
# Select data -------------------------------------------------------------####
  data_to_select <- c("Module", "Module_description", "Pathway", 
                      "Pathway_description", "Genes", 
                      "Gene_description", "Enzyme", "Cycle", "Pathway_cycle",
                      "Detail_cycle", "rbims_pathway", "rbims_sub_pathway", 
                      "KO", "dbCAN", "domain_name", "Pfam", "PFAM", "INTERPRO")
  
# Transform from wide to long ---------------------------------------------####
  Kegg_long <- tibble_ko %>%
    tidyr::pivot_longer(
      cols = tidyselect::where(is.numeric) | tidyselect::starts_with("Bin_"),
      values_to = "Abundance",
      names_to  = "Bin_name"
    ) %>%
    dplyr::distinct()
  # Count the total number of genes per metabolism---------------------------####
  metabolism_counts_totals <- Kegg_long %>%
    dplyr::select(!!y_axis_enquo, KO) %>%           # KO sin .data$
    dplyr::distinct() %>%
    dplyr::count(!!y_axis_enquo, sort = TRUE) %>%
    dplyr::rename(Abundance_metabolism = n) %>%     # n sin .data$
    dplyr::mutate(new2 = stringr::str_replace(Abundance_metabolism, "$", "\\)")) %>%
    dplyr::mutate(new2 = stringr::str_replace(new2, "^", "\\(")) %>%
    tidyr::unite("New", !!y_axis_enquo, new2, sep = " ", remove = FALSE) %>%  # nombres crudos
    dplyr::select(New, !!y_axis_enquo, Abundance_metabolism)
  
  # Count per genome --------------------------------------------------------####
  metabolism_counts_genome <- Kegg_long %>%
    dplyr::select(!!y_axis_enquo, Bin_name, KO, Abundance) %>%  # nombres crudos
    dplyr::distinct() %>%
    dplyr::filter(Abundance != 0) %>%                           # numérico, no "0"
    dplyr::group_by(Bin_name) %>%
    dplyr::count(!!y_axis_enquo) %>%
    dplyr::rename(Abundance_metabolism_genome = n)
  
# Calculate the percentage-------------------------------------------------####
  Table_with_percentage <- left_join(
    metabolism_counts_genome, metabolism_counts_totals, 
    by = y_axis_label) %>%
    mutate(Percentage = (.data$Abundance_metabolism_genome * 100) / 
                         .data$Abundance_metabolism ) %>%
    select(.data$New, {{y_axis_enquo}}, .data$Bin_name, .data$Percentage) %>%
    distinct() %>%
    drop_na() %>%
    rename(Pathway_number_of_total_elements = .data$New)
# Join data experiment ----------------------------------------------------####
  if(is.null(data_experiment) == F){
    Table_with_percentage<-Table_with_percentage %>%
    left_join(data_experiment, by = "Bin_name")
  }
  
  return(Table_with_percentage)
}