#' @title Metabolism table
#' @description Creates a metabolism table that includes the information from
#' KEGG and Interpro annotation.
#' @usage write_metabolism(data_interpro=NULL, data_ko=NULL)
#' @param data_interpro a table, output of InterProScan on tsv format.
#' InterProScan should have been run with -pa option to be able to use the 
#' KEGG option, in the database argument.
#' @param data_ko path to the directory where the KEGG output files are.
#' @details This function is part of a package used for the analysis of bins 
#' metabolism.
#' @import tibble dplyr stringr tidyr janitor rlang writexl
#' @importFrom purrr map reduce 
#' @examples
#' \dontrun{
#' write_metabolism("inst/extdata/Interpro_test.tsv", 
#' "inst/extdata/")
#' }
#' @export
write_metabolism<-function(
  data_interpro=NULL,
  data_ko=NULL){
  # Read the interpro -------------------------------------------------------####
  databases_names<-c("Pfam", "INTERPRO", "TIGRFAM", "SUPERFAMILY", "SMART", 
                     "SFLD", "ProSiteProfiles", "ProSitePatterns", "ProDom",
                     "PRINTS", "PIRSF", "MobiDBLite", "Hamap", "Gene3D",
                     "Coils", "CDD")
  All_interpro<-map(databases_names, read_interpro_colapsed, 
                    data_interpro=data_interpro) %>%
    reduce(full_join, by = c("Bin_name", "Scaffold_name")) %>%
    distinct()
  message("Finishing parsing Interpro")
  # Read KO data-------------------------------------------------------------####
  data_to_select<-c("Module", "Module_description", "Pathway", 
                    "Pathway_description", "Cycle", "Pathway_cycle",
                    "Detail_cycle", "Genes", "Gene_description", 
                    "Enzyme", "KO",
                    "rbims_pathway", "rbims_sub_pathway")
  ko_mapped<- read_ko(data_ko)  %>%
    mapping_ko() 
  ko<-read_ko(data_ko) %>%
    select(-.data$Abundance)  %>%
    group_by(.data$Bin_name, .data$Scaffold_name) %>%
    summarise(KO = paste0(.data$KO, collapse = ", "), 
              .groups = 'drop') %>%
    ungroup()
  All_data<-ko_mapped %>%
    pivot_longer(!all_of(data_to_select), names_to = "Bin_name", 
                 values_to = "Abundance") %>%
    select(-.data$Abundance) %>%
    distinct()
  message("Finishing parsing KEGG")
  # Creating main metabolism ----------------------------------------------####
  main_metabolism<-map(data_to_select, separate_kegg, 
                       kegg_data=ko, metabolism_data=All_data) %>%
    reduce(full_join, by = c("Bin_name", "Scaffold_name")) %>%
    full_join(All_interpro, by = c("Bin_name", "Scaffold_name")) %>%
    distinct()   %>%
    mutate(
      across(everything(), ~replace_na(.x, "---"))
    ) %>%
    mutate(
      across(everything(), ~str_replace_all(.x, "NA", "---"))
    ) %>%
    mutate(
      across(everything(), ~str_replace_all(.x, "---,", "---"))
    ) %>%
    mutate(
      across(everything(), ~str_replace_all(.x, "--- ---", "---"))
    ) %>%
    mutate(
      across(everything(), ~str_replace_all(.x, "--- --- ---", "---"))
    ) %>%
    mutate(
      across(everything(), ~str_replace_all(.x, "--- ---", "---"))
    )
  rm(All_interpro, ko)
  message("Finishing the Metabolism table")
  # Fixing mapping ko -----------------------------------------------------####
  ko_mapped<- ko_mapped %>%
    mutate(
      across(all_of(data_to_select), ~replace_na(.x, "---"))
    )
  ko<-read_ko(data_ko)
  # Read the interpro -------------------------------------------------------####
  message("Creating the Metabolism book")
  Pfam<-read_interpro(data_interpro, database="Pfam", profile = F) %>%
    rename(pfam_domain_name = .data$domain_name)   %>%
    mutate(Abundance=as.character(.data$Abundance))  %>%
    mutate(
      across(everything(), ~replace_na(.x, "---"))
    ) %>%
    mutate_at('Abundance', as.integer)
  INTERPRO<-read_interpro(data_interpro, database="INTERPRO", profile = F)%>%
    rename(INTERPRO_domain_name = .data$domain_name)  %>%
    mutate(Abundance=as.character(.data$Abundance))  %>%
    mutate(
      across(everything(), ~replace_na(.x, "---"))
    ) %>%
    mutate_at('Abundance', as.integer)
  TIGRFAM<-read_interpro(data_interpro, database="TIGRFAM", profile = F)%>%
    rename(TIGRFAM_domain_name = .data$domain_name)  %>%
    mutate(Abundance=as.character(.data$Abundance))  %>%
    mutate(
      across(everything(), ~replace_na(.x, "---"))
    ) %>%
    mutate_at('Abundance', as.integer)
  SUPERFAMILY<-read_interpro(data_interpro, database="SUPERFAMILY", 
                             profile = F) %>%
    rename(SUPERFAMILY_domain_name = .data$domain_name)  %>%
    mutate(Abundance=as.character(.data$Abundance))  %>%
    mutate(
      across(everything(), ~replace_na(.x, "---"))
    ) %>%
    mutate_at('Abundance', as.integer)
  SMART<-read_interpro(data_interpro, database="SMART", profile = F)%>%
    rename(SMART_domain_name = .data$domain_name)  %>%
    mutate(Abundance=as.character(.data$Abundance))  %>%
    mutate(
      across(everything(), ~replace_na(.x, "---"))
    ) %>%
    mutate_at('Abundance', as.integer)
  SFLD<-read_interpro(data_interpro, database="SFLD", profile = F)%>%
    rename(SFLD_domain_name = .data$domain_name)  %>%
    mutate(Abundance=as.character(.data$Abundance))  %>%
    mutate(
      across(everything(), ~replace_na(.x, "---"))
    ) %>%
    mutate_at('Abundance', as.integer)
  ProSiteProfiles<-read_interpro(data_interpro, database="ProSiteProfiles", 
                                 profile = F)%>%
    rename(ProSiteProfiles_domain_name = .data$domain_name)  %>%
    mutate(Abundance=as.character(.data$Abundance))  %>%
    mutate(
      across(everything(), ~replace_na(.x, "---"))
    ) %>%
    mutate_at('Abundance', as.integer)
  ProSitePatterns<-read_interpro(data_interpro, database="ProSitePatterns", 
                                 profile = F)%>%
    rename(ProSitePatterns_domain_name = .data$domain_name)  %>%
    mutate(Abundance=as.character(.data$Abundance))  %>%
    mutate(
      across(everything(), ~replace_na(.x, "---"))
    ) %>%
    mutate_at('Abundance', as.integer)
  ProDom<-read_interpro(data_interpro, database="ProDom", profile = F)%>%
    rename(ProDom_domain_name = .data$domain_name)   %>%
    mutate(Abundance=as.character(.data$Abundance))  %>%
    mutate(
      across(everything(), ~replace_na(.x, "---"))
    ) %>%
    mutate_at('Abundance', as.integer)
  PRINTS<-read_interpro(data_interpro, database="PRINTS", profile = F)%>%
    rename(PRINTS_domain_name = .data$domain_name)   %>%
    mutate(Abundance=as.character(.data$Abundance))  %>%
    mutate(
      across(everything(), ~replace_na(.x, "---"))
    ) %>%
    mutate_at('Abundance', as.integer)
  PIRSF<-read_interpro(data_interpro, database="PIRSF", profile = F)%>%
    rename(PIRSF_domain_name = .data$domain_name)  %>%
    mutate(Abundance=as.character(.data$Abundance))  %>%
    mutate(
      across(everything(), ~replace_na(.x, "---"))
    ) %>%
    mutate_at('Abundance', as.integer)
  MobiDBLite<-read_interpro(data_interpro, database="MobiDBLite", 
                            profile = F)%>%
    rename(MobiDBLite_domain_name = .data$domain_name)  %>%
    mutate(Abundance=as.character(.data$Abundance))  %>%
    mutate(
      across(everything(), ~replace_na(.x, "---"))
    ) %>%
    mutate_at('Abundance', as.integer)
  Hamap<-read_interpro(data_interpro, database="Hamap", profile = F)%>%
    rename(Hamap_domain_name = .data$domain_name)  %>%
    mutate(Abundance=as.character(.data$Abundance))  %>%
    mutate(
      across(everything(), ~replace_na(.x, "---"))
    ) %>%
    mutate_at('Abundance', as.integer)
  Gene3D<-read_interpro(data_interpro, database="Gene3D", profile = F)%>%
    rename(Gene3D_domain_name = .data$domain_name)  %>%
    mutate(Abundance=as.character(.data$Abundance))  %>%
    mutate(
      across(everything(), ~replace_na(.x, "---"))
    ) %>%
    mutate_at('Abundance', as.integer)
  Coils<-read_interpro(data_interpro, database="Coils", profile = F)%>%
    rename(Coils_domain_name = .data$domain_name)  %>%
    mutate(Abundance=as.character(.data$Abundance))  %>%
    mutate(
      across(everything(), ~replace_na(.x, "---"))
    ) %>%
    mutate_at('Abundance', as.integer)
  CDD<-read_interpro(data_interpro, database="CDD", profile = F)%>%
    rename(CDD_domain_name = .data$domain_name)  %>%
    mutate(Abundance=as.character(.data$Abundance))  %>%
    mutate(
      across(everything(), ~replace_na(.x, "---"))
    ) %>%
    mutate_at('Abundance', as.integer)
  # List sheets -------------------------------------------------------------####
  list_of_datasets <-
    list("MetabolismTable" = main_metabolism, 
         "KEGG profile" = ko_mapped, 
         "KEGG" = ko, 
         "Pfam" = Pfam,
         "InterProScanIDs" = INTERPRO,
         "TIGRFAM" = TIGRFAM,
         "SUPERFAMILY" = SUPERFAMILY,
         "SMART" = SMART,
         "SFLD" = SFLD,
         "ProSiteProfiles" = ProSiteProfiles, 
         "ProSitePatterns" = ProSitePatterns,
         "ProDom" = ProDom, 
         "PRINTS" = PRINTS, 
         "PIRSF" = PIRSF,
         "MobiDBLite" = MobiDBLite,
         "Hamap" = Hamap,
         "Gene3D" = Gene3D,
         "Coils" = Coils,
         "CDD" = CDD)
  # Write final table -------------------------------------------------------####
  write_xlsx(list_of_datasets, "MetabolismWorkbook.xlsx")
  message("Done")
}