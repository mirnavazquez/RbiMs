#' @title Read interpro databases.
#' @description read each interpro database.
#' @usage read_interpro_databases("Interprscan_tsv", database_1="Pfam")
#' @param data_interpro a table, output of InterProScan on tsv format.
#' InterProScan should have been run with -pa option to be able to use the 
#' KEGG option, in the database argument.
#' @param database_1 a character, indicating the name of the database 
#' to colapse.
#' @details This function is part of a package used for the analysis of bins 
#' metabolism.
#' @import  tibble dplyr stringr tidyr  rlang
#' @noRd
read_interpro_databases<- function(
  data_interpro,
  database_1){
  # Quote the database ----------------------------------------------------####
  database_quote<-enquo(database_1)
  # Make the table --------------------------------------------------------####
  Pfam<-read_interpro(data_interpro=data_interpro, database=database_1, 
                      profile = F) %>%
    mutate(
      across(everything(), ~replace_na(.x, "---"))
    ) %>%
    mutate_at('Abundance', as.integer)
  # Print -----------------------------------------------------------------####
  return(Pfam)
}