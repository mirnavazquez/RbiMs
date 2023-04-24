#' @title Subsets the mapped KO table based on pathways.
#' @description reads the output of the mapping_ko function to filter out KOs 
#' based on pathways.
#' @usage get_subset_pathway(mapped_ko_table, type_of_interest_feature,
#' interest_feature)
#' @param mapped_ko_table a data frame object. The output table from 
#' the mapping_ko function. 
#' @param type_of_interest_feature is a mapped_ko column name. This feature 
#' is going to be used for subsetting.  
#' @param interest_feature a character vector of the mapped_ko 
#' feature of interest. It is found under the type_of_interest_feature 
#' value column.
#' @details This function is part of a package used for the analysis of 
#' bins metabolism.
#' @import dplyr rlang
#' @examples
#' get_subset_pathway(mapped_ko_table=ko_bin_mapp, 
#' type_of_interest_feature=rbims_pathway, interest_feature="Hexadecane")
#' @export
get_subset_pathway<-function(mapped_ko_table,
                             type_of_interest_feature,
                             interest_feature){
  # Quoting ---------------------------------------------------------------####
  type_of_interest_feature_enquo <- enquo(type_of_interest_feature)
  # Filtering -------------------------------------------------------------####
  final_table<-mapped_ko_table %>%
    filter(!!type_of_interest_feature_enquo %in% interest_feature) %>%
    distinct()
  return(final_table)
}