library(dplyr)

get_unique_metabolism<-function(df_pfam, groupOfInterest){
  gruop_1<-filter(df_pfam, profundidad == groupOfInterest) %>%
    select(-profundidad)
  allGroupNames <- unique(df_pfam$profundidad)
  otherGroups <- allGroupNames[which(allGroupNames  != groupOfInterest )]
  The_rest <- filter(df_pfam, profundidad %in%  otherGroups) %>%
    select(-profundidad)
  uniquePFAMs <- c()
  for(i in 2:(length(gruop_1))){
    if(isTRUE(sum(gruop_1[,i]) > 0 & sum(The_rest[,i]) == 0)){
      uniquePFAMs <- c(uniquePFAMs, colnames(gruop_1[i]))
    }
  }
  return(uniquePFAMs)
}

b %>%
  select(-c(Module, Module_description, Pathway, 
         Pathway_description, Genes, 
         Gene_description, Enzyme, Scaffold_name)) %>%
  distinct() %>%
  column_to_rownames("KO") %>%
  as.data.frame(t(Final_master_sheet_short)) %>%  
  rownames_to_column("Bins") %>%
  left_join(metadata, by="Bin_name")

  
  
  
  
  
  
  
  
  
  
