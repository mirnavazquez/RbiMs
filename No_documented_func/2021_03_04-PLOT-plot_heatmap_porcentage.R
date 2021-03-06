library(pheatmap)

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
  #################### Transform from wide to long ####################
  Kegg_long<- tabble_ko %>%
    pivot_longer(cols = -c(Module, Module_description, Pathway, 
                           Pathway_description, Genes, 
                           Gene_description, Enzyme, KO), 
                 values_to = "Abundance",
                 names_to="Bin_name") %>%
    distinct()
  ######### Count the total number of genes per metabolism ########
  metabolism_counts_totals<-Kegg_long %>%
    select({{metadata_feature_enquo}}, KO) %>%
    distinct() %>%
    count(!!metadata_feature_enquo, sort=T) %>%
    rename(Abundance_metabolism=n)
  ## Count the total number of genes per metabolism in each genome ##
  metabolism_counts_genome<-Kegg_long %>%
    select({{metadata_feature_enquo}}, Bin_name, KO, Abundance) %>%
    distinct() %>%
    filter(Abundance != "0") %>%
    group_by(Bin_name) %>%
    count(!!metadata_feature_enquo) %>%
    rename(Abundance_metabolism_genome=n)
  ###################################################################
  tabble_ko_colnames<-colnames(tabble_ko)
  paths<-c("Module", "Module_description", "Pathway", 
           "Pathway_description", "Genes", 
           "Gene_description", "Enzyme", "KO")
  tabble_ko_colnames_2 <- tabble_ko_colnames[tabble_ko_colnames %in%
                                               paths]
  tabble_ko_colnames_3 <- tabble_ko_colnames_2[!tabble_ko_colnames_2 %in%
                                                 metadata_feature_label]
  ############### Join and calculate the percentage #################
  Table_with_percentage<- left_join(
    metabolism_counts_genome, metabolism_counts_totals, 
    by= metadata_feature_label) %>%
    mutate(Percentage = 
             (Abundance_metabolism_genome * 100) / 
             Abundance_metabolism ) %>%
    select({{metadata_feature_enquo}}, Bin_name, Percentage) %>%
    distinct() %>%
    drop_na() %>%
    pivot_wider(names_from = Bin_name, 
                values_from = Percentage, 
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
    column_to_rownames(metadata_feature_enquo)
  
  return(Table_with_percentage)
}

plot_heatmap_percentage(subset_mapped, metadata, Pathway, Modules,Clades)

    ################ Row metadata ##############
    metadata_hydro_path<- tabble_ko %>%
      select({{metadata_feature_enquo}}, {{row_feature_enquo}}) %>%
      drop_na() %>%
      distinct(Pathway_description, .keep_all = T) %>%
      column_to_rownames(metadata_feature_enquo) %>%
      arrange(!!col_feature_enquo)
    ################ Col metadata ##############
    metadata_column<- tabble_ko %>% 
      pivot_longer(cols = -c(Module, Module_description, Pathway, 
                             Pathway_description, Genes, 
                             Gene_description, Enzyme, KO), 
                   values_to = "Abundance",
                   names_to="Bin_name") %>%
      distinct() %>%
      left_join(other_data, by="Bin_name") %>%
      select( Bin_name, {{row_feature_enquo}}) %>%
      drop_na() %>%
      distinct(Bin_name, .keep_all = T) %>%
      column_to_rownames("Bin_name") %>%
      arrange(!!row_feature_enquo)
    
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
    
  
  
  
  