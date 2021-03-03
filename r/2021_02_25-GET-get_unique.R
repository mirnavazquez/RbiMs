#' @title Identifies unique KOs in the data
#' @description reads the output of the mapping_ko function 
#' to filter out KOs that are not mapped to a unique protein family.
#' @param mapped_ko_table is the output table from the mapping_ko function. 
#' @param other_data is a table that contains the necessary metadata. 
#' @param type_of_interest_feature is the metadata column name containing different features. 
#' @param interest_feature is a string of the metadata feature of interest found under the 
#' type_of_interest_feature value column.
#' @details This function is part of a package used for the analysis of bins metabolism.
#' @import dplyr
#' @examples
#' get_unique(kegg_mapped, metadata, Sample_site, "Water_column")
#' @export
get_unique<-function(mapped_ko_table, 
                     other_data,
                     type_of_interest_feature, 
                     interest_feature){
  ############################ quoting ######################################
  all_features <- enquo(type_of_interest_feature)
  all_features_x <- as_label(all_features)
  ############################# Parse the mapped ko table ###################
  mapped_ko_table_KO<-mapped_ko_table %>%
    select(-c(Module, Module_description, Pathway, 
              Pathway_description, Genes, 
              Gene_description, Enzyme)) %>%
    distinct() %>%
    column_to_rownames("KO") 
  ############################# Transpose the table ########################
  mapped_ko_table_KO_t<-as.data.frame(t(mapped_ko_table_KO)) %>%  
    rownames_to_column("Bin_name") %>%
    left_join(other_data, by="Bin_name")
  ############################# Get metadata names ##########################
  metadata_names<-colnames(other_data)
  metadata_names<-metadata_names[-1]
  ############################# Run the unique lines  #######################
  gruop_1<-mapped_ko_table_KO_t%>%
    filter( !!all_features  == interest_feature) %>%
    select(-any_of(metadata_names))
  allGroupNames <- mapped_ko_table_KO_t %>%
    select(all_of(all_features_x)) %>%
    distinct()
  otherGroups <- allGroupNames %>%
    filter(!!all_features != interest_feature)
  otherGroups_list<-lapply(otherGroups, as.character)
  otherGroups_nolist<-unlist(otherGroups_list)
  The_rest <- mapped_ko_table_KO_t %>%
    filter( !!all_features %in% otherGroups_nolist) %>%
    select(-all_of(metadata_names))
  ########################## loop to get the uniques  #######################
  uniquePFAMs <- c()
  for(i in 2:(length(gruop_1))){
    if(isTRUE(sum(gruop_1[,i]) > 0 & sum(The_rest[,i]) == 0)){
      uniquePFAMs <- c(uniquePFAMs, colnames(gruop_1[i]))
    }
  }
  ########################## Extract the values #############################
  final_table <- mapped_ko_table %>%
    filter(KO %in% uniquePFAMs)
  return(final_table)
}