#' @title Calculate presence/absence.
#' @description  Calculate the  presence/absence of KO in certain pathway.
#' @usage calc_binary(tibble_ko, y_axis, data_experiment=NULL, binary=TRUE,
#' metabolism=FALSE)
#' @param tibble_ko a tibble object from mapping_ko.
#' @param y_axis a character, indicating the pathway to analyze.
#' @param data_experiment optional. a data frame object 
#' containing metadata information.
#' @param binary optional. a logical value, indicating whether you want to
#' calculate presence/absence or abundance of KOs per pathway.
#' @param metabolism optional. a logical value, if TRUE it will show 
#' the metabolism table.
#' @details This function is part of a package used for the analysis of bins 
#' metabolism.
#' @import  tibble dplyr stringr tidyr  rlang
#' @examples
#' calc_binary(ko_bin_mapp, Pathway)    
#' @export
calc_binary<-function(tibble_ko,
                      y_axis,
                      data_experiment=NULL,
                      binary=TRUE,
                      metabolism=FALSE){
  # Enquoting -------------------------------------------------------------####
  y_axis_enquo <- enquo(y_axis)
  y_axis_label <- as_label(y_axis_enquo)
  # Select data -----------------------------------------------------------####
  data_to_select<-c("Module", "Module_description", "Pathway", 
                    "Pathway_description", "Genes", 
                    "Gene_description", "Enzyme", "Cycle", "Pathway_cycle",
                    "Detail_cycle", "rbims_pathway", "rbims_sub_pathway",
                    "KO")
  # Transform from wide to long -------------------------------------------####
  Kegg_long<- tibble_ko %>%
    pivot_longer(cols = -all_of(data_to_select), 
                 values_to = "Abundance",
                 names_to="Bin_name") %>%
    distinct()
  # Calc binary -----------------------------------------------------------####
  if(binary == TRUE){
    Table_binary<-suppressMessages(Kegg_long %>%
                                     rename(tmp = .data$Abundance)  %>%
                                     select({{y_axis_enquo}}, .data$Bin_name, 
                                            .data$tmp) %>%
                                     distinct() %>%
                                     group_by( {{y_axis_enquo}}, 
                                               .data$Bin_name)  %>%
                                     summarise(Presence_absence = 
                                                 sum(.data$tmp)) %>% 
                                     ungroup() %>%
                                     drop_na() %>%
                                     mutate(Presence_absence = case_when(
                                       .data$Presence_absence != "0" ~ "1",
                                       TRUE ~ as.character(.data$Presence_absence)
                                     )) %>%
                                     mutate_at('Presence_absence', as.integer)) 
  } else if (binary == FALSE) {
    Table_binary<-suppressMessages(Kegg_long  %>%
                                     select({{y_axis_enquo}}, .data$Bin_name, 
                                            .data$Abundance) %>%
                                     distinct() %>%
                                     group_by({{y_axis_enquo}}, 
                                              .data$Bin_name) %>%
                                     summarise(Abundance = 
                                                 sum(.data$Abundance)) %>% 
                                     ungroup() %>%
                                     drop_na())
  }
  # Join data experiment --------------------------------------------------####
  if(is.null(data_experiment) == F){
    Table_binary<-Kegg_long %>%
      left_join(Table_binary, by="Bin_name")
  }
  # Full table ----------------------------------------------------------- ####
  if(metabolism == T){
    Kegg_long<-select(Kegg_long, -.data$Abundance)
    Table_binary<-Table_binary %>%
      left_join(Kegg_long, by=c(y_axis_label, "Bin_name"))
  }
  return(Table_binary)
}