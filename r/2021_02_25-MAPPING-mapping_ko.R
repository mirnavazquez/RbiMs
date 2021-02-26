#' @title mapping_ko 
#' @description It reads an object created with read_ko and maps each KO
#' to the KEGG database.
#' @param ko_abundance_table an object created with read_ko.
#' @details This function is part of a package used for 
#' the analysis of bins metabolism.
#' @import KEGGREST tibble dplyr stringr tidyr janitor
#' @examples
#' mapping_ko(kegg_read)
#' @export
mapping_ko<-function(ko_abundance_table){
  #########################  KEGG links ############################
  module_link<- keggLink("ko", "module") %>%
    enframe() %>%
    rename(KO = value) %>%
    rename(Module = name) %>%
    mutate(Module=str_remove_all(Module, "md:")) %>%
    mutate(KO=str_remove_all(KO, "ko:")) 
  pathway_link <- keggLink("ko", "pathway") %>%
    enframe() %>%
    rename(KO = value) %>%
    rename(Pathway = name) %>%
    filter(!str_detect(Pathway, "ko")) %>%
    mutate(Pathway=str_remove_all(Pathway, "path:")) %>%
    mutate(KO=str_remove_all(KO, "ko:"))  
  enzyme_link<-keggLink("ko", "enzyme") %>%
    enframe() %>%
    rename(KO = value) %>%
    rename(Enzyme = name ) %>%
    mutate(KO=str_remove_all(KO, "ko:")) 
  #########################  KEGG lists ############################
  modules_list<-keggList("module") %>%
    enframe() %>%
    rename(Module_description = value) %>%
    rename(Module = name ) %>%
    mutate(Module=str_remove_all(Module, "md:")) 
  pathway_list<-keggList("pathway")%>%
    enframe() %>%
    rename(Pathway_description = value) %>%
    rename(Pathway = name ) %>%
    mutate(Pathway=str_remove_all(Pathway, "path:")) 
  ko_list<-keggList("ko") %>%
    enframe() %>%
    rename(KO_description = value) %>%
    rename(KO = name ) %>%
    separate(KO_description, 
           c("Genes", "Gene_description"), sep="; ") %>%
    mutate(KO=str_remove_all(KO, "ko:"))
  #########################  KO master #############################
  KO_master<-left_join(ko_list, pathway_link, by="KO") %>%
    left_join(pathway_list, by="Pathway") %>%
    left_join(module_link, by="KO") %>%
    left_join(modules_list, by="Module") %>%
    left_join(enzyme_link, by="KO") %>%
    select(Module, Module_description, Pathway, 
         Pathway_description, KO, Genes, 
         Gene_description, Enzyme) %>%
    distinct()
  ######################  Remove data ##############################
  rm(module_link, pathway_link, enzyme_link, modules_list, 
    pathway_list, ko_list)
  ###################################################################
  ko_abundance_table %>%
  left_join(KO_master, by="KO") %>%
    select(KO, Bin_name, Abundance) %>%
    distinct() %>%
  left_join(KO_master, by="KO") %>%
    select(Module, Module_description, Pathway, 
         Pathway_description, KO, Abundance, Genes, 
         Gene_description, Enzyme,Bin_name) %>%
    distinct() %>%
    pivot_wider(names_from = Bin_name, 
                values_from = Abundance, values_fill = 0) 
}





