#' @title Separate KEGG to add to the database.
#' @description Divide the KEGG database.
#' @usage separate_kegg(kegg_data=kegg_data, 
#' metabolism_data=metabolism_data, database=database)
#' @param kegg_data a table with KEGG information.
#' @param metabolism_data a table created within write_metabolism.
#' @param database a character. The feature of the KEGG db to parse.
#' @details This function is part of a package used for the analysis of bins 
#' metabolism.
#' @import  tibble dplyr stringr tidyr  rlang
#' @noRd
separate_kegg<- function(
  kegg_data,
  metabolism_data,
  database){
  # Quote the database ----------------------------------------------------####
  database_quote<-enquo(database)
  # Make the table --------------------------------------------------------####
  Module<-kegg_data %>%
    left_join(metabolism_data, by =c("Bin_name", "KO"))  %>%
    group_by(.data$Bin_name, .data$Scaffold_name) %>%
    select(.data$Bin_name, .data$Scaffold_name, {{database_quote}}) %>%
    distinct() %>%
    summarise(!!database_quote := paste0(.data[[database]], collapse = ", "), 
              .groups = 'drop') %>%
    ungroup()
  # Print -----------------------------------------------------------------####
  return(Module)
}