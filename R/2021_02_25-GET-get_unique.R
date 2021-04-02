#' @title Mapped KO table sub-setting. Finding the KOs only present 
#' in a specific metadata group.
#' @description reads the output of the mapping_ko function to pull out KOs 
#' based on metadata. This function identifies the KOs that are only present 
#' in a particular metadata group (i.e. KOs only present in soil and 
#' absent in the rest of the environments). And retrieves a tibble object which 
#' is a sub-set of the mapped_ko table.
#' @param tibble_ko_mapp a tibble object. The output table from 
#' the mapping_ko function. 
#' @param data_experiment a data frame object that contains the metadata 
#' (i.e. taxonomy, sampling site). 
#' @param experiment_col a metadata column name. 
#' This feature is going to be used for sub-setting.  
#' @param experiment_col_element a string of the metadata feature of interest. 
#' It is found under the experiment_col column.
#' @details This function is part of a package used for the analysis 
#' of bins metabolism.
#' @import dplyr rlang
#' @examples
#' get_subset_unique(ko_bin_mapp, metadata, Sample_site, "Water_column")
#' @export
get_subset_unique<-function(tibble_ko_mapp, 
                            data_experiment,
                            experiment_col, 
                            experiment_col_element){
  # Variable quoting ------------------------------------------------------####
  experiment_col_enquo <- enquo(experiment_col)
  experiment_col_label <- as_label(experiment_col_enquo)
  # Remove columns, create df_unique ------------------------------------####
  data_to_select<-c("Module", "Module_description", "Pathway", 
                    "Pathway_description", "Genes", 
                    "Gene_description", "Enzyme", "Cycle", "Pathway_cycle",
                    "Detail_cycle",  "rbims_pathway", "rbims_sub_pathway")
  df_unique<-tibble_ko_mapp %>%
    select(-all_of(data_to_select)) %>%
    distinct() %>%
    column_to_rownames("KO")
  # Columns to variables, merge with experiment ---------------------------####
  df_unique_t<-as.data.frame(t(df_unique)) %>%  
    rownames_to_column("Bin_name") %>%
    left_join(data_experiment, by="Bin_name")
  # Extract experiment names ----------------------------------------------####
  experiment_names<-colnames(data_experiment)
  experiment_names<-experiment_names[-1]
  # Extract bins based on experiment_col_element --------------------------####
  df_exp_col_element<-df_unique_t%>%
    filter( !!experiment_col_enquo  == experiment_col_element) %>%
    select(-any_of(experiment_names))
  # Extract bins based not in experiment_col_element ----------------------####
  vector_all_exp_col_val <- df_unique_t %>%
    select(all_of(experiment_col_label)) %>%
    distinct() %>%
    filter(!!experiment_col_enquo != experiment_col_element) %>%
    pull()
  df_no_exp_col_element <- df_unique_t %>%
    filter( !!experiment_col_enquo %in% vector_all_exp_col_val) %>%
    select(-all_of(experiment_names))
  # Loop data -------------------------------------------------------------####
  list_unique <- c()
  for(i in 2:(length(df_exp_col_element))){
    if(isTRUE(sum(df_exp_col_element[,i]) > 0 & sum(df_no_exp_col_element[,i]) == 0)){
      list_unique <- c(list_unique, colnames(df_exp_col_element[i]))
    }
  }
  # Write tibble ----------------------------------------------------------####
  final_table <- tibble_ko_mapp %>%
    filter(.data$KO %in% list_unique)
  return(final_table)
}