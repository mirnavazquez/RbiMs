#' @title Maps KOs to the KEGG database
#' @description Reads a table object created with the read_ko function 
#' and maps each KO to the KEGG database.
#' @param ko_abundance_table a data frame object with four columns containing the KO abundance 
#' for each Bin.
#' @details This function is part of a package used for the analysis of bins metabolism.
#' @import KEGGREST tibble dplyr stringr tidyr janitor rlang
#' @examples
#' mapping_ko(ko_bin_table)
#' @export
mapping_ko<-function(ko_abundance_table){
  #########################  KEGG links ############################
  module_link<- KEGGREST::keggLink("ko", "module") %>%
    enframe() %>%
    rename(KO = .data$value) %>%
    rename(Module = .data$name) %>%
    mutate(Module=str_remove_all(.data$Module, "md:")) %>%
    mutate(KO=str_remove_all(.data$KO, "ko:")) 
  pathway_link <- keggLink("ko", "pathway") %>%
    enframe() %>%
    rename(KO = .data$value) %>%
    rename(Pathway = .data$name) %>%
    filter(!str_detect(.data$Pathway, "ko")) %>%
    mutate(Pathway=str_remove_all(.data$Pathway, "path:")) %>%
    mutate(KO=str_remove_all(.data$KO, "ko:"))  
  enzyme_link<-keggLink("ko", "enzyme") %>%
    enframe() %>%
    rename(KO = .data$value) %>%
    rename(Enzyme = .data$name ) %>%
    mutate(KO=str_remove_all(.data$KO, "ko:")) 
  #########################  KEGG lists ############################
  modules_list<-KEGGREST::keggList("module") %>%
    enframe() %>%
    rename(Module_description = .data$value) %>%
    rename(Module = .data$name ) %>%
    mutate(Module=str_remove_all(.data$Module, "md:")) 
  pathway_list<-KEGGREST::keggList("pathway")%>%
    enframe() %>%
    rename(Pathway_description = .data$value) %>%
    rename(Pathway = .data$name ) %>%
    mutate(Pathway=str_remove_all(.data$Pathway, "path:")) 
  ko_list<-suppressWarnings(keggList("ko") %>%
                              enframe() %>%
                              rename(KO_description = .data$value) %>%
                              rename(KO = .data$name ) %>%
                              separate(.data$KO_description, 
                                       c("Genes", "Gene_description"), sep="; ") %>%
                              mutate(KO=str_remove_all(.data$KO, "ko:")))
  #########################  KO master #############################
  data_to_select<-c("Module", "Module_description", "Pathway", 
                    "Pathway_description", "KO", "Genes", "Gene_description", "Enzyme") 
  KO_master<-left_join(ko_list, pathway_link, by="KO") %>%
    left_join(pathway_list, by="Pathway") %>%
    left_join(module_link, by="KO") %>%
    left_join(modules_list, by="Module") %>%
    left_join(enzyme_link, by="KO") %>%
    select(all_of(data_to_select)) %>%
    distinct()
  ######################  Remove data ##############################
  rm(module_link, pathway_link, enzyme_link, modules_list, 
     pathway_list, ko_list, data_to_select)
  ######################  Add DiTing ###############################
  DiTing_cycles<-suppressMessages(read_delim(
    "https://raw.githubusercontent.com/xuechunxu/DiTing/master/table/KO_affilated_to_biogeochemical_cycle.tab", 
    delim="\t") %>%
      fill(.data$Cycle) %>%
      fill(.data$Pathway) %>%
      rename(Pathway_cycle = .data$Pathway) %>%
      rename(KO = .data$k_number) %>%
      rename(Detail_cycle = .data$Detail))
  #####################  Combine datasets ###########################
  KO_master_DiTing<-left_join(KO_master, DiTing_cycles, by="KO")
  ###################################################################
  data_to_select<-c("Module", "Module_description", "Pathway", 
                    "Pathway_description", "Cycle", "Pathway_cycle",
                    "Detail_cycle", "Genes", "Gene_description", 
                    "Enzyme", "Bin_name", "Abundance")
  final_table <- ko_abundance_table %>%
    left_join(KO_master_DiTing, by="KO") %>%
    select(.data$KO, .data$Bin_name, .data$Abundance) %>%
    distinct() %>%
    left_join(KO_master_DiTing, by="KO") %>%
    select(all_of(data_to_select)) %>%
    distinct() %>%
    pivot_wider(names_from = .data$Bin_name, 
                values_from = .data$Abundance, values_fill = 0) 
  return(final_table)
}