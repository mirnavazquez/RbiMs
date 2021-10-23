#' @title Colapse the interpro databases.
#' @description Colapse the interpro databases.
#' @usage read_interpro_colapsed(data_interpro=data_interpro, 
#' database_1=database_1)
#' @param data_interpro a table, output of InterProScan on tsv format.
#' InterProScan should have been run with -pa option to be able to use the 
#' KEGG option, in the database argument.
#' @param database_1 a character, indicating the name of the database 
#' to colapse.
#' @details This function is part of a package used for the analysis of bins 
#' metabolism.
#' @import  tibble dplyr stringr tidyr  rlang
#' @noRd
read_interpro_colapsed<- function(
  data_interpro,
  database_1){
  # Quote the database ----------------------------------------------------####
  database_quote<-enquo(database_1)
  # Make the table --------------------------------------------------------####
  Pfam<-read_interpro(data_interpro=data_interpro, database=database_1, 
                      profile = F) %>%
    rename(pfam_domain_name = .data$domain_name)  %>%
    group_by(.data$Bin_name, .data$Scaffold_name) %>%
    select(.data$Bin_name, .data$Scaffold_name, {{database_quote}})   %>%
    distinct() %>%
    summarise(!!database_quote := paste0(.data[[database_1]], collapse = ", "), 
              .groups = 'drop') %>%
    ungroup()
  # Print -----------------------------------------------------------------####
  return(Pfam)
}