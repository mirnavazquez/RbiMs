#' @title Heatmap plot of KO abundance within each Bin.
#' @description Creates a heatmap of KO abundance.
#' @param tabble_ko a data frame object, created with the mapping_ko or get_subset_* functions.
#' @param other_data a data frame object containing metadata information.
#' @param plot_ano is a column name indicating the annotation if rows. This column name comes from the mapping_ko or get_subset_* object.
#' @param plot_medata is a column name indicating the annotation if cols. This column name comes from the metadata object.
#' @details This function is part of a package used for the analysis of bins metabolism.
#' @import pheatmap rlang
#' @examples
#' plot_heatmap_abundance(ko_bin_mapp, metadata, Module, Sample_site)
#' @export
plot_heatmap_abundance<-function(tabble_ko, 
                                 other_data,
                                 plot_ano,
                                 plot_medata){
  ############################ quoting ##############################
  plot_ano_enquo <- enquo(plot_ano)
  plot_medata_enquo <- enquo(plot_medata)
  plot_ano_label <- as_label(plot_ano_enquo)
  ###################################################################
  tabble_ko_colnames<-colnames(tabble_ko)
  paths<-c("Module", "Module_description", "Pathway", 
           "Pathway_description", "Genes", 
           "Gene_description", "Enzyme", "KO", "Cycle", "Pathway_cycle",
           "Detail_cycle")
  tabble_ko_colnames_2 <- tabble_ko_colnames[tabble_ko_colnames %in% paths]
  tabble_ko_colnames_3 <- tabble_ko_colnames_2[!tabble_ko_colnames_2 %in% "KO"]
  both_paths<-c("KO", plot_ano_label)
  tabble_ko_colnames_4 <- tabble_ko_colnames_3[!tabble_ko_colnames_3 %in% plot_ano_label]
  ########## Create pivot tables #############
  metabolism.wide <- tabble_ko %>%
    select(-all_of(tabble_ko_colnames_3)) %>%
    distinct() %>%
    column_to_rownames("KO")
  ################ Row metadata ##############
  metadata_hydro_path<- tabble_ko %>%
    select( "KO", {{ plot_ano_enquo }}) %>%
    drop_na() %>%
    distinct(.data$KO, .keep_all = T) %>%
    column_to_rownames("KO") %>%
    arrange({{ plot_ano_enquo }})
  ################ Row metadata ##############
  data_to_select<-c("Module", "Module_description", "Pathway", 
                    "Pathway_description", "Genes", 
                    "Gene_description", "Enzyme", "Cycle", "Pathway_cycle",
                    "Detail_cycle", "KO")
  metadata_column<- tabble_ko %>% 
    pivot_longer(cols = -data_to_select, values_to = "Abundance",
                 names_to="Bin_name") %>%
    distinct() %>%
    left_join(other_data, by="Bin_name") %>%
    select( "Bin_name", {{ plot_medata_enquo }}) %>%
    drop_na() %>%
    distinct(.data$Bin_name, .keep_all = T) %>%
    column_to_rownames("Bin_name") %>%
    arrange({{ plot_medata_enquo }})
  #######  Reorder the names of rows ##############################
  sub_samp_ordered <- metabolism.wide[rownames(metadata_hydro_path),]
  sub_samp_ordered_2 <- sub_samp_ordered[,rownames(metadata_column)]
  
  metabolism.wide_pheatmap<-pheatmap(sub_samp_ordered_2, 
                                     scale = "row",
                                     annotation_row = metadata_hydro_path,
                                     annotation_col = metadata_column,
                                     cluster_rows = F,
                                     cluster_cols = F)
  
  return(metabolism.wide_pheatmap)
}