library(pheatmap)


tabble_ko_colnames<-colnames(subset_mapped)
paths<-c("Module", "Module_description", "Pathway", 
         "Pathway_description", "Genes", 
         "Gene_description", "Enzyme", "KO")
tabble_ko_colnames_2 <- tabble_ko_colnames[tabble_ko_colnames %in% paths]
tabble_ko_colnames_3 <- tabble_ko_colnames_2[!tabble_ko_colnames_2 %in% "Pathway"]
both_paths<-c("KO", plot_ano_label)
tabble_ko_colnames_4 <- tabble_ko_colnames_3[!tabble_ko_colnames_3 %in% plot_ano_label]

  #################### Transform from wide to long ####################
  Kegg_long<- subset_mapped %>%
    pivot_longer(cols = -c(Module, Module_description, Pathway, 
                           Pathway_description, Genes, 
                           Gene_description, Enzyme, KO), 
                 values_to = "Abundance",
                 names_to="Bin_name") %>%
    distinct()
  ######### Count the total number of genes per metabolism ########
  metabolism_counts_totals<-Kegg_long %>%
    select(Pathway, KO) %>%
    distinct() %>%
    count(Pathway, sort=T) %>%
    rename(Abundance_metabolism=n)
  ## Count the total number of genes per metabolism in each genome ##
  metabolism_counts_genome<-Kegg_long %>%
    select(Pathway, Bin_name, KO, Abundance) %>%
    distinct() %>%
    filter(Abundance != "0") %>%
    group_by(Bin_name) %>%
    count(Pathway) %>%
    rename(Abundance_metabolism_genome=n)
  ############### Join and calculate the percentage #################
  Table_with_percentage<- left_join(
    metabolism_counts_genome, metabolism_counts_totals, 
    by= "Pathway") %>%
    mutate(Percentage = 
             (Abundance_metabolism_genome * 100) / Abundance_metabolism ) %>%
    select(Pathway, Bin_name, Percentage) %>%
    distinct() %>%
    drop_na() %>%
    pivot_wider(names_from = Bin_name, 
                values_from = Percentage, 
                values_fill = 0) %>%
    distinct() %>%
    left_join(subset_mapped, 
              by= "Pathway") 
  
  str_match(Table_with_percentage, ".y")  
  
  a<-Table_with_percentage %>%
  select(-c(Module, Module_description,
            Pathway_description, Genes, 
              Gene_description, Enzyme, KO), -contains(".y")) %>%
    distinct() %>%
    rename_all(
      list( ~ stringr::str_replace_all(., ".x", ""))
    ) %>%
    column_to_rownames("Pathway")
  
  
  return(Table_with_percentage)
}


plot_heatmap_percentage(subset_mapped, metadata, Pathway, Modules,Clades)
################ Row metadata ##############
metadata_hydro_path<- tabble_ko %>%
  select(Pathway, {{row_feature_enquo}}) %>%
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




