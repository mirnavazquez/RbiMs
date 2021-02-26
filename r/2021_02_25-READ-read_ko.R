#' @title read_ko read the output of KofamScan or KofamKoala
#' @description read the output of KofamScan or KofamKoala and 
#' calculate the abundance of each KO within the Bins.
#' @param ko_result KEGG output file. Contings are expected 
#' to be named similar to bin_scaffoldXX.
#' @details This function is part of a package used for 
#' the analysis of bins metabolism.
#' @import dplyr tidyr readr stringr
#' @examples
#' read_ko("KEGG_bins.txt")
#' @export
read_ko<-function(ko_result){
  KO_raw<-read_table2(ko_result, col_names = F) %>%
    filter(!str_detect(X1, '#')) %>%
    select(X2, X3) %>%
    rename(Bin_name = X2) %>%
    rename(KO = X3) %>%
    separate(Bin_name, c("Bin_name", "Scaffold_name"), 
             sep = "[_|-][s|S]caffold") %>%
    mutate(Scaffold_name = paste0( "scaffold", Scaffold_name), 
           Scaffold_name) 
  
  KO_abundance<-KO_raw %>%
    group_by(Bin_name) %>%
    distinct() %>%
    count(KO) %>%
    rename(Abundance = n)%>%
    ungroup()
  
  left_join(KO_raw, KO_abundance, by=c("Bin_name", "KO")) %>%
    distinct()

}
  



