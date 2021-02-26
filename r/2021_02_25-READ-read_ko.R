library(dplyr)
library(tidyr)
library(readr)
library(stringr)

read_ko<-function(ko_result){
  KO_raw<-read_table2(ko_result, col_names = F) %>%
    filter(!str_detect(X1, '#')) %>%
    select(X2, X3) %>%
    rename(Bin_name = X2) %>%
    rename(KO = X3) %>%
    separate(Bin_name, c("Bin_name", "Scaffold_name"), sep = "[_|-][s|S]caffold") %>%
    mutate(Scaffold_name = paste0( "scaffold", Scaffold_name), Scaffold_name) 
  
  KO_abundance<-KO_raw %>%
    group_by(Bin_name) %>%
    distinct() %>%
    count(KO) %>%
    rename(Abundance = n)%>%
    ungroup()
  
  left_join(KO_raw, KO_abundance, by=c("Bin_name", "KO")) %>%
    distinct()

}
  
kegg_read<-read_ko("data/KEGG_bins.txt")


