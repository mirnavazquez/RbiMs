#' @title Heatmap plot of KO/Pathways/Modules and its relative 
#' abundance within each bin.
#' @description Creates a heatmap of KO/Pathways/Modules and 
#' its relative abundance within each Bin.
#' @param tibble_ko a data frame object, created with the mapping_ko 
#' or get_subset_* functions.
#' @param data_experiment a data frame object containing metadata information.
#' @param metabolism_col_sort a string. A column name of the tibble_ko of a 
#' feature to plot (i.e. KO/Pathways/Modules).
#' @param metabolism_col is a column name indicating the annotation if rows. 
#' This column name comes from the mapping_ko or get_subset_* object.
#' @param experiment_col is a column name indicating the annotation if cols. 
#' This column name comes from the metadata object.
#' @details This function is part of a package used for 
#' the analysis of bins metabolism.
#' @import pheatmap rlang dplyr tidyr tibble
#' @examples
#' plot_heatmap_percentage(ko_bin_mapp, metadata, Pathway, Module, Clades) 
#' @export
plot_heatmap_percentage<-function(tibble_ko,
                                  data_experiment, 
                                  metabolism_col_sort,
                                  metabolism_col,
                                  experiment_col){
  # Quoting ---------------------------------------------------------------####
  metabolism_col_sort_enquo <- enquo(metabolism_col_sort)
  metabolism_col_enquo <- enquo(metabolism_col)
  experiment_col_enquo <- enquo(experiment_col)
  metabolism_col_sort_label <- as_label(metabolism_col_sort_enquo)
  metabolism_col_label <- as_label(metabolism_col_enquo)
  # Transform from wide to long -------------------------------------------####
  data_to_select<-c("Module", "Module_description", "Pathway", 
                    "Pathway_description", "Genes", 
                    "Gene_description", "Enzyme", "KO", "Cycle", 
                    "Pathway_cycle",
                    "Detail_cycle", "rbims_pathway", "rbims_sub_pathway")
  tibble_ko_long<- tibble_ko %>%
    pivot_longer(cols = -data_to_select, 
                 values_to = "Abundance",
                 names_to="Bin_name") %>%
    distinct()
  # Count the total number of genes per metabolism ------------------------####
  metabolism_counts_totals<-tibble_ko_long %>%
    select({{metabolism_col_sort_enquo}}, .data$KO) %>%
    distinct() %>%
    count(!!metabolism_col_sort_enquo, sort=T) %>%
    rename(Abundance_metabolism=n)
  # Count the total number of genes per metabolism in each genome ---------####
  metabolism_counts_genome<-tibble_ko_long %>%
    select({{metabolism_col_sort_enquo}}, .data$Bin_name, .data$KO, 
           .data$Abundance) %>%
    distinct() %>%
    filter(.data$Abundance != "0") %>%
    group_by(.data$Bin_name) %>%
    count(!!metabolism_col_sort_enquo) %>%
    rename(Abundance_metabolism_genome=n)
  # Sorting col names -----------------------------------------------------####
  tibble_ko_colnames<-colnames(tibble_ko)
  paths<-c("Module", "Module_description", "Pathway", 
           "Pathway_description", "Genes", 
           "Gene_description", "Enzyme", "KO", "Cycle", "Pathway_cycle",
           "Detail_cycle", "rbims_pathway", "rbims_sub_pathway")
  tibble_ko_colnames_2 <- tibble_ko_colnames[tibble_ko_colnames %in%
                                               paths]
  tibble_ko_colnames_3 <- tibble_ko_colnames_2[!tibble_ko_colnames_2 %in%
                                                 metabolism_col_sort_label]
  # Join and calculate the percentage -------------------------------------####
  Table_with_percentage<- left_join(
    metabolism_counts_genome, metabolism_counts_totals, 
    by= metabolism_col_sort_label) %>%
    mutate(Percentage = 
             (.data$Abundance_metabolism_genome * 100) / 
             .data$Abundance_metabolism ) %>%
    select({{metabolism_col_sort_enquo}}, .data$Bin_name, 
           .data$Percentage) %>%
    distinct() %>%
    drop_na() %>%
    pivot_wider(names_from = .data$Bin_name, 
                values_from = .data$Percentage, 
                values_fill = 0) %>%
    distinct() %>%
    left_join(tibble_ko, 
              by= metabolism_col_sort_label) %>%
    select(-all_of(tibble_ko_colnames_3)) %>%
    select(-contains(".y")) %>%
    distinct() %>%
    rename_all(
      list( ~ stringr::str_replace_all(., ".x", ""))
    ) %>%
    column_to_rownames(metabolism_col_sort_label)
  # Extracting metabolism -------------------------------------------------####
  metabolism_order<- tibble_ko %>%
    select({{metabolism_col_sort_enquo}}, {{metabolism_col_enquo}}) %>%
    drop_na() %>%
    distinct(.data[[metabolism_col_sort_enquo]], .keep_all = T) %>%
    column_to_rownames(metabolism_col_sort_label) %>%
    arrange(!!metabolism_col_enquo)
  # Extracting experiment -------------------------------------------------####
  data_to_select<-c("Module", "Module_description", "Pathway", 
                    "Pathway_description", "Genes", 
                    "Gene_description", "Enzyme", "KO", "Cycle", 
                    "Pathway_cycle",
                    "Detail_cycle", "rbims_pathway", "rbims_sub_pathway")
  experiment_order<- tibble_ko %>% 
    pivot_longer(cols = -data_to_select, 
                 values_to = "Abundance",
                 names_to="Bin_name") %>%
    distinct() %>%
    left_join(data_experiment, by="Bin_name") %>%
    select( .data$Bin_name, {{experiment_col_enquo}}) %>%
    drop_na() %>%
    distinct(.data$Bin_name, .keep_all = T) %>%
    column_to_rownames("Bin_name") %>%
    arrange(!!experiment_col_enquo)
  # Order table -----------------------------------------------------------####
  sub_samp_ordered <- Table_with_percentage[rownames(metabolism_order),]
  sub_samp_ordered_2 <- sub_samp_ordered[,rownames(experiment_order)]
  # Plot ------------------------------------------------------------------####
  metabolism.wide_pheatmap<-pheatmap(sub_samp_ordered_2, 
                                     scale = "none",
                                     annotation_row = metabolism_order,
                                     annotation_col = experiment_order,
                                     cluster_rows = F,
                                     cluster_cols = F)
  
  return(metabolism.wide_pheatmap)
}  
  
  
  