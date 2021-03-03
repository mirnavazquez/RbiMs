#' @title Obtains PCA data of KO pathways 
#' @description Identifies the most important KO pathways in each Bin.
#' @param ko_table is a table object created with the get_unique function
#' and contains unique KOs.
#' @details This function is part of a package used for 
#' the analysis of bins metabolism.
#' @import dplyr factoextra
#' @examples
#' get_pca_metabolism(Kegg_mapping)
#' @export
get_pca_metabolism<-function(ko_table){
##################### Read the input table ########################
  wide_ko<-ko_table %>%
    select(-c(Module, Module_description, Pathway, 
              Pathway_description, Genes, 
              Gene_description, Enzyme)) %>%
    distinct() %>%
    column_to_rownames("KO") 
  ################## Calculate the distance #######################
  wider_dist<-dist(wide_ko)
  #################### Calculate the distance #####################
  df_pca <- prcomp(wider_dist, center = T, scale = T)
  #################### Extract PCA infromation ####################
  pca_information<-get_pca(df_pca)
  #################### Extract PCA infromation ####################
  contribution_Metabolism<-as.data.frame(pca_information$cos2)
  ##### Warning if the contribution <=0.98 ########################
  if(contribution_Metabolism$Dim.1 <= 0.98) 
    warning("Contribution of Dim.1 is less than 0.98. Your output
            will be the same as input")
  ##### Create subsets of the first components ####################
  subset1_pathways<-subset(contribution_Metabolism, 
                         Dim.1 >= 0.98)
  subset2_pathways<-subset(contribution_Metabolism, 
                         Dim.2  >= 0.98)
  subset_098<-rbind(subset1_pathways, 
                  subset2_pathways)
  KO_important_paths<- subset(wide_ko, 
                      rownames(wide_ko) %in% rownames(subset_098))
  KO_important_paths_2<-rownames_to_column(Kegg_mapping,
                                         var = "ko")
  return(KO_important_paths_2)
}







