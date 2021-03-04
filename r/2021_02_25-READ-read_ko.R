#' @title Reads the output of KofamScan or KofamKoala
#' @description Calculates the abundance of each KO within the 
#' Bins of KofamScan or KofamKoala data.
#' @param ko_result is a string of the KEGG output table file. 
#' Contigs are expected to be named similar to bin_scaffoldXX.
#' @details This function is part of a package used for the analysis of bins metabolism.
#' @import dplyr tidyr readr stringr
#' @examples
#' read_ko("KEGG_bins.txt")
#' @export
read_ko<-function(ko_result){
  KO_raw<-suppressWarnings(read_table2(ko_result, col_names = F) %>%
    filter(!str_detect(X1, '#')) %>%
    select(X2, X3) %>%
    rename(Bin_name = X2) %>%
    rename(KO = X3) %>%
    separate(Bin_name, c("Bin_name", "Scaffold_name"), 
             sep = "[_|-][s|S]caffold") %>%
    mutate(Scaffold_name = paste0( "scaffold", Scaffold_name), 
           Scaffold_name))
  
  KO_abundance<-KO_raw %>%
    group_by(Bin_name) %>%
    distinct() %>%
    count(KO) %>%
    rename(Abundance = n)%>%
    ungroup()
  
  final_table <- left_join(KO_raw, KO_abundance, by=c("Bin_name", "KO")) %>%
    distinct()
  
  return(final_table)
}
  



