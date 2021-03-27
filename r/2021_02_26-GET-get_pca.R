#' @title Metabolism PCA. 
#' @description Identifies the most important KO pathways in the whole database.
#' @param tibble_ko_map a tibble object, created with 
#' the mapping_ko or get_subset_* functions. 
#' @param cos2_val a numeric vector from 0 to 1 indicating the proportion of
#' contribution used as cut off. Default is 0.98. 
#' See \link[factoextra]{get_pca}
#' @details This function is part of a package used for 
#' the analysis of bins metabolism.
#' @import dplyr factoextra rlang tibble
#' @examples
#' get_subset_pca(ko_bin_mapp)
#' @export
get_subset_pca<-function(tibble_ko_map,
                         cos2_val=NULL){
  # Read data -------------------------------------------------------------####
  data_to_select<-c("Module", "Module_description", "Pathway", 
                    "Pathway_description", "Genes", 
                    "Gene_description", "Enzyme", "Cycle", "Pathway_cycle",
                    "Detail_cycle", "rbims_pathway", "rbims_sub_pathway")
  wide_ko<-tibble_ko_map %>%
    dplyr::select(-all_of(data_to_select)) %>%
    dplyr::distinct() %>%
    tibble::column_to_rownames("KO") 
  # Distance --------------------------------------------------------------####
  wider_dist<-stats::dist(wide_ko)
  # PCA -------------------------------------------------------------------####
  df_pca <- stats::prcomp(wider_dist, center = T, scale = T)
  # Extract PCA contribution ----------------------------------------------####
  pca_information<-factoextra::get_pca(df_pca)
  contribution_Metabolism<-as.data.frame(pca_information$cos2)
  # Warning if the contribution <=0.98 ------------------------------------####
  if (all(contribution_Metabolism$Dim.1 <= 0.97)){
    warning("Contribution of the first dimention is less or equal to 0.97")
  }else { 
    message("Contribution of the first dimention is higher than 0.97")
  } 
  # Extract PC1,2 ---------------------------------------------------------####
  if(is.null(cos2_val) == T){
    cos2_val <- 0.98
  }
  # Extract PC1,2 ---------------------------------------------------------####
  subset1_pathways<-subset(contribution_Metabolism, 
                           contribution_Metabolism$Dim.1 >= cos2_val)
  subset2_pathways<-subset(contribution_Metabolism, 
                           contribution_Metabolism$Dim.2 >= cos2_val)
  subset1<-c(rownames(subset1_pathways), rownames(subset2_pathways))
  # Write tibble ----------------------------------------------------------####
  Final_table<-tibble_ko_map %>%
    filter(.data$KO %in% subset1)
  
  return(Final_table)
}






