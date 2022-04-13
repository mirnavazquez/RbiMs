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
                         analysis=c("KEGG", "Pfam", "INTERPRO")){
  # Read table ------------------------------------------------------------####
  KO_raw<-tabla_toabundance %>%
    separate(.data$Bin_name, c("Bin_name", "Scaffold_name"),
             sep = "[_|-][s|S]caffold") %>%
    mutate(Scaffold_name = paste0( "scaffold", .data$Scaffold_name),
           .data$Scaffold_name) %>%
    unite("Scaffold_name", c("Bin_name", "Scaffold_name"), remove=FALSE)
  # Selecting analysis ----------------------------------------------------####
  if( analysis == "Pfam") {
    KO_raw <- KO_raw %>% 
      rename(tmp = .data$Pfam) %>%
      distinct()
  } else if ( analysis == "KEGG") {
    KO_raw <- KO_raw %>% 
      rename(tmp = .data$KO)
  } else if ( analysis == "INTERPRO"){
    KO_raw <- KO_raw %>% 
      rename(tmp = .data$Interpro) %>%
      distinct()
  } else if ( analysis == "TIGRFAM"){
    KO_raw <- KO_raw %>% 
      rename(tmp = .data$TIGRFAM) %>%
      distinct()
  } else if ( analysis == "SUPERFAMILY"){
    KO_raw <- KO_raw %>% 
      rename(tmp = .data$SUPERFAMILY) %>%
      distinct()
  } else if ( analysis == "SMART"){
    KO_raw <- KO_raw %>% 
      rename(tmp = .data$SMART) %>%
      distinct()
  } else if ( analysis == "SFLD"){
    KO_raw <- KO_raw %>% 
      rename(tmp = .data$SFLD) %>%
      distinct()
  } else if ( analysis == "ProSiteProfiles"){
    KO_raw <- KO_raw %>% 
      rename(tmp = .data$ProSiteProfiles) %>%
      distinct()
  } else if ( analysis == "ProSitePatterns"){
    KO_raw <- KO_raw %>% 
      rename(tmp = .data$ProSitePatterns) %>%
      distinct()
  } else if ( analysis == "ProDom"){
    KO_raw <- KO_raw %>% 
      rename(tmp = .data$ProDom) %>%
      distinct()
  } else if ( analysis == "PRINTS"){
    KO_raw <- KO_raw %>% 
      rename(tmp = .data$PRINTS) %>%
      distinct()
  } else if ( analysis == "PIRSF"){
    KO_raw <- KO_raw %>% 
      rename(tmp = .data$PIRSF) %>%
      distinct()
  } else if ( analysis == "MobiDBLite"){
    KO_raw <- KO_raw %>% 
      rename(tmp = .data$MobiDBLite) %>%
      distinct()
  } else if ( analysis == "Hamap"){
    KO_raw <- KO_raw %>% 
      rename(tmp = .data$Hamap) %>%
      distinct()
  } else if ( analysis == "Gene3D"){
    KO_raw <- KO_raw %>% 
      rename(tmp = .data$Gene3D) %>%
      distinct()
  } else if ( analysis == "Coils"){
    KO_raw <- KO_raw %>% 
      rename(tmp = .data$Coils) %>%
      distinct()
  } else if ( analysis == "CDD"){
    KO_raw <- KO_raw %>% 
      rename(tmp = .data$CDD) %>%
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
  
  if( analysis == "Pfam") {
    final_table <- final_table_1 %>% 
      rename(Pfam = .data$tmp) 
  } else if ( analysis == "KEGG")  {
    final_table <- final_table_1 %>% 
      rename(KO = .data$tmp)
  } else if ( analysis == "INTERPRO"){
    final_table <- final_table_1 %>% 
      rename(INTERPRO = .data$tmp)
  } else if ( analysis == "TIGRFAM"){
    final_table <- final_table_1 %>% 
      rename(TIGRFAM = .data$tmp)
  } else if ( analysis == "SUPERFAMILY"){
    final_table <- final_table_1 %>% 
      rename(SUPERFAMILY = .data$tmp)
  } else if ( analysis == "SMART"){
    final_table <- final_table_1 %>% 
      rename(SMART = .data$tmp)
  } else if ( analysis == "SFLD"){
    final_table <- final_table_1 %>% 
      rename(SFLD = .data$tmp)
  } else if ( analysis == "ProSiteProfiles"){
    final_table <- final_table_1 %>% 
      rename(ProSiteProfiles = .data$tmp)
  } else if ( analysis == "ProSitePatterns"){
    final_table <- final_table_1 %>% 
      rename(ProSitePatterns = .data$tmp)
  } else if ( analysis == "ProDom"){
    final_table <- final_table_1 %>% 
      rename(ProDom = .data$tmp)
  } else if ( analysis == "PRINTS"){
    final_table <- final_table_1 %>% 
      rename(PRINTS = .data$tmp)
  } else if ( analysis == "PIRSF"){
    final_table <- final_table_1 %>% 
      rename(PIRSF = .data$tmp)
  } else if ( analysis == "MobiDBLite"){
    final_table <- final_table_1 %>% 
      rename(MobiDBLite = .data$tmp)
  } else if ( analysis == "Hamap"){
    final_table <- final_table_1 %>% 
      rename(Hamap = .data$tmp)
  } else if ( analysis == "Gene3D"){
    final_table <- final_table_1 %>% 
      rename(Gene3D = .data$tmp)
  } else if ( analysis == "Coils"){
    final_table <- final_table_1 %>% 
      rename(Coils = .data$tmp)
  } else if ( analysis == "CDD"){
    final_table <- final_table_1 %>% 
      rename(CDD = .data$tmp)
  } 
  
  return(final_table)
} 