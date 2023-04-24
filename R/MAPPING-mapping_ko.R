#' @title Maps KOs to the KEGG database.
#' @description Reads a table object created with the read_ko function 
#' and maps each KO to the KEGG database.
#' @usage mapping_ko(tibble_ko=NULL, tibble_interpro=NULL)
#' @param tibble_ko a data frame object with four columns containing 
#' the KO abundance for each bin. Output of read_ko.
#' @param tibble_interpro a data frame with three columns. The first column  
#' refers to contigs. They are expected to indicate in their names
#' the bin name follow by the scaffold name divided by an "_": 
#' bin_scaffoldXX. The second column is the KEGG pathway to which
#' such scaffold belong, and the third is the enzymes. 
#' Output of read_interpro.
#' @details This function is part of a package used for the analysis of bins 
#' metabolism.
#' @import KEGGREST tibble dplyr stringr tidyr janitor rlang
#' @examples
#' mapping_ko(ko_bin_table)
#' @export
mapping_ko<-function(tibble_ko=NULL, 
                     tibble_interpro=NULL ){
  # KEGG links ------------------------------------------------------------####
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
  # KEGG lists ------------------------------------------------------------####
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
  ko_list<-suppressMessages(suppressWarnings(
    read_delim("http://rest.kegg.jp/list/ko",  
               delim="\t", col_names = F) %>%
    rename(KO_description = .data$X2) %>%
    rename(KO = .data$X1 ) %>%
    separate(.data$KO_description, 
             c("Genes", "Gene_description"), 
             sep="; ",  extra = "merge") %>%
    mutate(KO=str_remove_all(.data$KO, "ko:")))) 
  
  #ko_list<-suppressWarnings(keggList("ko") %>%
   #                           enframe() %>%
    #                          rename(KO_description = .data$value) %>%
     #                         rename(KO = .data$name ) %>%
      #                        separate(.data$KO_description, 
       #                                c("Genes", "Gene_description"), 
        #                               sep="; ") %>%
         #                     mutate(KO=str_remove_all(.data$KO, "ko:")))
  # KO master -------------------------------------------------------------####
  data_to_select<-c("Module", "Module_description", "Pathway", 
                    "Pathway_description", "KO", "Genes", 
                    "Gene_description", "Enzyme") 
  KO_master<-left_join(ko_list, pathway_link, by="KO") %>%
    left_join(pathway_list, by="Pathway") %>%
    left_join(module_link, by="KO") %>%
    left_join(modules_list, by="Module") %>%
    left_join(enzyme_link, by="KO") %>%
    select(all_of(data_to_select)) %>%
    distinct()
  # Remove data -----------------------------------------------------------####
  rm(module_link, pathway_link, enzyme_link, modules_list, 
     pathway_list, ko_list, data_to_select)
  # Add DiTing ------------------------------------------------------------####
  DiTing_cycles<-suppressMessages(read_delim(
    "https://raw.githubusercontent.com/xuechunxu/DiTing/master/table/KO_affilated_to_biogeochemical_cycle.tab", 
    delim="\t") %>%
      fill(.data$Cycle) %>%
      fill(.data$Pathway) %>%
      rename(Pathway_cycle = .data$Pathway) %>%
      rename(KO = .data$k_number) %>%
      rename(Detail_cycle = .data$Detail))
  # Combine data-sets -----------------------------------------------------####
  KO_master_DiTing<-left_join(KO_master, DiTing_cycles, by="KO")
  # Read RbiMs ------------------------------------------------------------####
  rbims<-rbims
  # Combine data-sets -----------------------------------------------------####
  KO_master_DiTing_rbims<-left_join(KO_master_DiTing, rbims, by="KO") %>%
    mutate(Genes = case_when(
      Genes == "alkT" ~ "rubB, alkT",
      TRUE ~ as.character(Genes)
    ))
  # Data to select --------------------------------------------------------####
  data_to_select<-c("Module", "Module_description", "Pathway", 
                    "Pathway_description", "Cycle", "Pathway_cycle",
                    "Detail_cycle", "Genes", "Gene_description", 
                    "Enzyme", "Bin_name", "Abundance", "KO",
                    "rbims_pathway", "rbims_sub_pathway")
  
  if(is.null(tibble_interpro) == F){
    output_interpro <- tibble_interpro %>%
      left_join(KO_master_DiTing_rbims, by=c("Pathway", "Enzyme"))
    Interpro_ko<-read_ko(data_kofam=NULL, data_kaas=NULL, output_interpro)
    table_to_ko<- Interpro_ko
  }
  if(is.null(tibble_ko) == F){
    table_to_ko<- tibble_ko
  }
  # Write tibble ----------------------------------------------------------####
  final_table <- table_to_ko %>%
    left_join(KO_master_DiTing_rbims, by="KO") %>%
    select(.data$KO, .data$Bin_name, .data$Abundance) %>%
    distinct() %>%
    left_join(KO_master_DiTing_rbims, by="KO") %>%
    select(all_of(data_to_select)) %>%
    distinct() %>%
    pivot_wider(names_from = .data$Bin_name, 
                values_from = .data$Abundance, values_fill = 0)
  
  return(final_table)
}
