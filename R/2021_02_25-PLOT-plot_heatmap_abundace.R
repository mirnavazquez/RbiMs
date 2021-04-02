#' @title Heatmap plot of KO abundance within each bin.
#' @description Creates a heatmap of KO abundance.
#' @param tibble_ko a data frame object, created with the mapping_ko or 
#' get_subset_* functions.
#' @param data_experiment a data frame object containing metadata information.
#' @param metabolism_col is a column name indicating the annotation if rows. 
#' This column name comes from the mapping_ko or get_subset_* object.
#' @param experiment_col is a column name indicating the annotation of cols. 
#' This column name comes from the metadata object.
#' @details This function is part of a package used for 
#' the analysis of bins metabolism.
#' @import pheatmap rlang
#' @examples
#' plot_heatmap_abundance(ko_bin_mapp, metadata, Module, Sample_site)
#' @export
plot_heatmap_abundance<-function(tibble_ko, 
                                 data_experiment,
                                 metabolism_col,
                                 experiment_col){
  # Quoting ---------------------------------------------------------------####
  metabolism_col_enquo <- enquo(metabolism_col)
  experiment_col_enquo <- enquo(experiment_col)
  metabolism_col_label <- as_label(metabolism_col_enquo)
  # Sorting col names -----------------------------------------------------####
  tibble_ko_colnames<-colnames(tibble_ko)
  paths<-c("Module", "Module_description", "Pathway", 
           "Pathway_description", "Genes", 
           "Gene_description", "Enzyme", "KO", "Cycle", "Pathway_cycle",
           "Detail_cycle", "rbims_pathway", "rbims_sub_pathway")
  tibble_ko_colnames_2 <- tibble_ko_colnames[tibble_ko_colnames %in% paths]
  tibble_ko_colnames_3 <- tibble_ko_colnames_2[!tibble_ko_colnames_2 %in% "KO"]
  both_paths<-c("KO", metabolism_col_label)
  tibble_ko_colnames_4 <- tibble_ko_colnames_3[!tibble_ko_colnames_3 
                                               %in% metabolism_col_label]
  # Removing cols ---------------------------------------------------------####
  heatmap_table <- tibble_ko %>%
    select(-all_of(tibble_ko_colnames_3)) %>%
    distinct() %>%
    column_to_rownames("KO")
  # Extracting metabolism -------------------------------------------------####
  metabolism_order<- tibble_ko %>%
    select( "KO", {{ metabolism_col_enquo }}) %>%
    drop_na() %>%
    distinct(.data$KO, .keep_all = T) %>%
    column_to_rownames("KO") %>%
    arrange({{ metabolism_col_enquo }})
  # Extracting experiment -------------------------------------------------####
  data_to_select<-c("Module", "Module_description", "Pathway", 
                    "Pathway_description", "Genes", 
                    "Gene_description", "Enzyme", "Cycle", "Pathway_cycle",
                    "Detail_cycle", "KO", "rbims_pathway", "rbims_sub_pathway")
  experiment_order<- tibble_ko %>% 
    pivot_longer(cols = -data_to_select, values_to = "Abundance",
                 names_to="Bin_name") %>%
    distinct() %>%
    left_join(data_experiment, by="Bin_name") %>%
    select( "Bin_name", {{ experiment_col_enquo }}) %>%
    drop_na() %>%
    distinct(.data$Bin_name, .keep_all = T) %>%
    column_to_rownames("Bin_name") %>%
    arrange({{ experiment_col_enquo }})
  # Order table -----------------------------------------------------------####
  sub_samp_ordered <- heatmap_table[rownames(metabolism_order),]
  sub_samp_ordered_2 <- sub_samp_ordered[,rownames(experiment_order)]
  # Plot ------------------------------------------------------------------####
  heatmap_table_pheatmap<-pheatmap(sub_samp_ordered_2, 
                                     scale = "row",
                                     annotation_row = metabolism_order,
                                     annotation_col = experiment_order,
                                     cluster_rows = F,
                                     cluster_cols = F)
  
  return(heatmap_table_pheatmap)
}
