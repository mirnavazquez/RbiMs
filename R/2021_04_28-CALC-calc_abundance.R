#' @title Calculate abundance.
#' @description Calculate the abundance of certain entry, based on the
#' number of times it appears.
#' @usage calc_abundance(tabla_toabundance, 
#' analysis=c("KEGG", "PFAM", "INTERPRO"))
#' @param tabla_toabundance a data frame object with the values to be 
#' calculated.
#' @param analysis a character, indicating which analysis it is doing.
#' @details This function is part of a package used for the analysis of bins 
#' metabolism.
#' @import  tibble dplyr stringr tidyr  rlang
#' @noRd
calc_abundance<-function(tabla_toabundance,
                         analysis=c("KEGG", "PFAM", "INTERPRO")){
  # Read table ------------------------------------------------------------####
  KO_raw<-tabla_toabundance %>%
    separate(.data$Bin_name, c("Bin_name", "Scaffold_name"),
             sep = "[_|-][s|S]caffold|-S") %>%
    mutate(Scaffold_name = paste0( "scaffold", .data$Scaffold_name),
           .data$Scaffold_name) 
  # Selecting analysis ----------------------------------------------------####
  if( analysis == "PFAM") {
    KO_raw <- KO_raw %>% 
      rename(tmp = .data$pfam) %>%
      #select(-.data$Scaffold_name)%>%
      distinct()
  } else if ( analysis == "KEGG") {
    KO_raw <- KO_raw %>% 
      rename(tmp = .data$KO)
  } else if ( analysis == "INTERPRO"){
    KO_raw <- KO_raw %>% 
      rename(tmp = .data$Interpro) %>%
      #select(-.data$Scaffold_name)%>%
      distinct()
  }
  # Calculate abundance ---------------------------------------------------####
  KO_abundance<-KO_raw %>%
    group_by(.data$Bin_name) %>%
    distinct() %>%
    count(.data$tmp) %>%
    rename(Abundance = .data$n)%>%
    ungroup()
  # Write tibble -----------------------------------------------------------####
  final_table_1 <- left_join(KO_raw, KO_abundance, 
                             by=c("Bin_name", "tmp")) %>%
    distinct() 
  
  if( analysis == "PFAM") {
    final_table <- final_table_1 %>% 
      rename(PFAM = .data$tmp) 
  } else if ( analysis == "KEGG")  {
    final_table <- final_table_1 %>% 
      rename(KO = .data$tmp)
  } else if ( analysis == "INTERPRO"){
    final_table <- final_table_1 %>% 
      rename(INTERPRO = .data$tmp)
  }
  
  return(final_table)
} 