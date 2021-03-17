#' @title Heatmap plot of KO/Pathways/Modules and its relative abundance within each Bin.
#' @description Creates a heatmap of KO/Pathways/Modules and its relative abundance within each Bin.
#' @param tabble_ko a data frame object, created with the mapping_ko or get_subset_* functions.
#' @param other_data a data frame object containing metadata information.
#' @param metadata_feature a string. A column name of the tabble_ko of a feature to plot (i.e. KO/Pathways/Modules).
#' @param row_feature is a column name indicating the annotation if rows. This column name comes from the mapping_ko or get_subset_* object.
#' @param col_feature is a column name indicating the annotation if cols. This column name comes from the metadata object.
#' @details This function is part of a package used for the analysis of bins metabolism.
#' @import pheatmap rlang dplyr tidyr tibble
#' @examples
#' plot_heatmap_percentage(ko_bin_mapp, metadata, Pathway, Module, Clades) 
#' @export
plot_heatmap_percentage<-function(tabble_ko,
                                  other_data, 
                                  metadata_feature,
                                  row_feature,
                                  col_feature){
  ############################ quoting ##############################
  metadata_feature_enquo <- enquo(metadata_feature)
  row_feature_enquo <- enquo(row_feature)
  col_feature_enquo <- enquo(col_feature)
  metadata_feature_label <- as_label(metadata_feature_enquo)
  row_feature_label <- as_label(row_feature_enquo)
  #################### Transform from wide to long ####################
  data_to_select<-c("Module", "Module_description", "Pathway", 
                    "Pathway_description", "Genes", 
                    "Gene_description", "Enzyme", "KO", "Cycle", "Pathway_cycle",
                    "Detail_cycle", "rbims_pathway", "rbims_sub_pathway")
  Kegg_long<- tabble_ko %>%
    pivot_longer(cols = -data_to_select, 
                 values_to = "Abundance",
                 names_to="Bin_name") %>%
    distinct()
  ######### Count the total number of genes per metabolism ########
  metabolism_counts_totals<-Kegg_long %>%
    select({{metadata_feature_enquo}}, .data$KO) %>%
    distinct() %>%
    count(!!metadata_feature_enquo, sort=T) %>%
    rename(Abundance_metabolism=n)
  ## Count the total number of genes per metabolism in each genome ##
  metabolism_counts_genome<-Kegg_long %>%
    select({{metadata_feature_enquo}}, .data$Bin_name, .data$KO, .data$Abundance) %>%
    distinct() %>%
    filter(.data$Abundance != "0") %>%
    group_by(.data$Bin_name) %>%
    count(!!metadata_feature_enquo) %>%
    rename(Abundance_metabolism_genome=n)
  ###################################################################
  tabble_ko_colnames<-colnames(tabble_ko)
  paths<-c("Module", "Module_description", "Pathway", 
           "Pathway_description", "Genes", 
           "Gene_description", "Enzyme", "KO", "Cycle", "Pathway_cycle",
           "Detail_cycle", "rbims_pathway", "rbims_sub_pathway")
  tabble_ko_colnames_2 <- tabble_ko_colnames[tabble_ko_colnames %in%
                                               paths]
  tabble_ko_colnames_3 <- tabble_ko_colnames_2[!tabble_ko_colnames_2 %in%
                                                 metadata_feature_label]
  ############### Join and calculate the percentage #################
  Table_with_percentage<- left_join(
    metabolism_counts_genome, metabolism_counts_totals, 
    by= metadata_feature_label) %>%
    mutate(Percentage = 
             (.data$Abundance_metabolism_genome * 100) / 
             .data$Abundance_metabolism ) %>%
    select({{metadata_feature_enquo}}, .data$Bin_name, .data$Percentage) %>%
    distinct() %>%
    drop_na() %>%
    pivot_wider(names_from = .data$Bin_name, 
                values_from = .data$Percentage, 
                values_fill = 0) %>%
    distinct() %>%
    left_join(tabble_ko, 
              by= metadata_feature_label) %>%
    select(-all_of(tabble_ko_colnames_3)) %>%
    select(-contains(".y")) %>%
    distinct() %>%
    rename_all(
      list( ~ stringr::str_replace_all(., ".x", ""))
    ) %>%
    column_to_rownames(metadata_feature_label)
  ################ Row metadata ##############
  metadata_hydro_path<- tabble_ko %>%
    select({{metadata_feature_enquo}}, {{row_feature_enquo}}) %>%
    drop_na() %>%
    distinct(.data[[metadata_feature_enquo]], .keep_all = T) %>%
    column_to_rownames(metadata_feature_label) %>%
    arrange(!!row_feature_enquo)
  ################ Col metadata ##############
  data_to_select<-c("Module", "Module_description", "Pathway", 
                    "Pathway_description", "Genes", 
                    "Gene_description", "Enzyme", "KO", "Cycle", "Pathway_cycle",
                    "Detail_cycle", "rbims_pathway", "rbims_sub_pathway")
  metadata_column<- tabble_ko %>% 
    pivot_longer(cols = -data_to_select, 
                 values_to = "Abundance",
                 names_to="Bin_name") %>%
    distinct() %>%
    left_join(other_data, by="Bin_name") %>%
    select( .data$Bin_name, {{col_feature_enquo}}) %>%
    drop_na() %>%
    distinct(.data$Bin_name, .keep_all = T) %>%
    column_to_rownames("Bin_name") %>%
    arrange(!!col_feature_enquo)
  #######  Reorder the names of rows and cols ##############################
  sub_samp_ordered <- Table_with_percentage[rownames(metadata_hydro_path),]
  sub_samp_ordered_2 <- sub_samp_ordered[,rownames(metadata_column)]
  
  metabolism.wide_pheatmap<-pheatmap(sub_samp_ordered_2, 
                                     scale = "none",
                                     annotation_row = metadata_hydro_path,
                                     annotation_col = metadata_column,
                                     cluster_rows = F,
                                     cluster_cols = F)
  
  return(metabolism.wide_pheatmap)
}  
  
  
  