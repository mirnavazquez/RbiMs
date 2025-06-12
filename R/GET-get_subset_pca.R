#' @title Metabolism PCA. 
#' @description Identifies the most important KO pathways or protein 
#' domains in the whole database. And print back a profile of the protein 
#' domains that have higher contributions.
#' @usage get_subset_pca(tibble_rbims, cos2_val=NULL, 
#' analysis=c("KEGG", "Pfam", "INTERPRO"))
#' @param tibble_rbims a tibble object, created with 
#' the read_interpro, mapping_ko or get_subset_* functions. 
#' @param cos2_val a numeric vector from 0 to 1 indicating the proportion of
#' contribution used as cut off. Default is 0.98.
#' See \link[factoextra]{get_pca}.
#' @param analysis a character, indicating from which input do you want to
#' get the abundance profile. Valid options are "KEGG", "Pfam" or "INTERPRO". 
#' @details This function is part of a package used for
#' the analysis of bins metabolism.
#' @import dplyr factoextra rlang tibble
#' @examples
#' #get_subset_pca(ko_bin_mapp, analysis="KEGG")
#' @export
get_subset_pca<-function(tibble_rbims,
                         cos2_val=NULL,
                         analysis=c("KEGG", "Pfam", "INTERPRO", "dbCAN", 
                                    "MEROPS")) {
                         
                      
  if(analysis == "PFAM"){
    stop("P letter must be in capital followed by lower-case letters e.g 'Pfam'")
   }
  
  # Select data -------------------------------------------------------####
  if(analysis=="KEGG") {
    data_to_select<-c("Module", "Module_description", "Pathway", 
                      "Pathway_description", "Genes", 
                      "Gene_description", "Enzyme", "Cycle", "Pathway_cycle",
                      "Detail_cycle", "rbims_pathway", "rbims_sub_pathway")
  } else if (analysis=="Pfam" || analysis=="INTERPRO" || analysis=="dbCAN") {
    data_to_select<-"domain_name"
  }
  # Rename ----------------------------------------------------------------####
  if( analysis == "Pfam") {
    tibble_rbims <- tibble_rbims %>% 
      rename(tmp = .data$Pfam) 
  } else if ( analysis == "KEGG") {
    tibble_rbims <- tibble_rbims %>% 
      rename(tmp = .data$KO)
  } else if ( analysis == "INTERPRO"){
    tibble_rbims <- tibble_rbims %>% 
      rename(tmp = .data$INTERPRO)
  } else if ( analysis == "dbCAN"){
    tibble_rbims <- tibble_rbims %>% 
      rename(tmp = .data$dbCAN_family)
 }
  # Read data -------------------------------------------------------------####
  wide_ko<-tibble_rbims %>%
    dplyr::select(-all_of(data_to_select)) %>%
    dplyr::distinct() %>%
    tibble::column_to_rownames("tmp") 

  # PCA -------------------------------------------------------------------####
  if (!requireNamespace("FactoMineR", quietly = TRUE)) {
  install.packages("FactoMineR")
}
  df_pca <- FactoMineR::PCA(wide_ko, scale.unit = TRUE, graph = FALSE)
  contribution_Metabolism<-as.data.frame(df_pca$ind$cos2)
  
  # Warning if the contribution <=0.98 ------------------------------------####
  if (all(contribution_Metabolism$Dim.1 <= 0.97)){
    warning("Contribution of the first dimention is less or equal to 0.97.
            if the data frame has 0 observations, try a different cos2 value.")
  }else { 
    message("Contribution of the first dimention is higher than 0.97.
            if the data frame has 0 observations, try a different cos value.")
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
  final_table_1<-tibble_rbims %>%
    filter(.data$tmp %in% subset1)
  
  if( analysis == "Pfam") {
    final_table <- final_table_1 %>% 
      rename(Pfam = .data$tmp) 
  } else if ( analysis == "KEGG")  {
    final_table <- final_table_1 %>% 
      rename(KO = .data$tmp)
  } else if ( analysis == "INTERPRO"){
    final_table <- final_table_1 %>% 
      rename(INTERPRO = .data$tmp)
  } else if ( analysis == "dbCAN"){
    final_table <- final_table_1 %>% 
      rename(dbCAN_fam = .data$tmp)
  }
  return(final_table)
}


