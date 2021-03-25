#' @title Read the output of KofamScan or KofamKoala
#' @description read_ko() calculates the abundance of each KO within the 
#' bins based on the KofamScan or KofamKoala output.
#' @param data_ko a data frame with 5 columns. The contigs are expected to 
#' indicate in their names the bin name and then the scaffold name divided 
#' by an "_", similar to bin_scaffoldXX.
#' @details This function is part of a package used for the analysis 
#' of bins metabolism.
#' @import dplyr tidyr readr stringr rlang
#' @examples
#' \dontrun{
#' read_ko("KEGG_bins.txt")
#' }
#' @export
read_ko<-function(data_ko){
  # Read data -------------------------------------------------------------####
  KO_raw<-suppressWarnings(suppressMessages(read_table2(data_ko, 
                                                        col_names = F) %>%
    filter(!str_detect(.data$X1, '#')) %>%
    select(.data$X2, .data$X3) %>%
    rename(Bin_name = .data$X2) %>%
    rename(KO = .data$X3) %>%
    separate(.data$Bin_name, c("Bin_name", "Scaffold_name"), 
             sep = "[_|-][s|S]caffold|-S") %>%
    mutate(Scaffold_name = paste0( "scaffold", .data$Scaffold_name), 
           .data$Scaffold_name)))
  # Calculate abundance ---------------------------------------------------####
  KO_abundance<-KO_raw %>%
    group_by(.data$Bin_name) %>%
    distinct() %>%
    count(.data$KO) %>%
    rename(Abundance = .data$n)%>%
    ungroup()
  # Write tibble -----------------------------------------------------------####
  final_table <- left_join(KO_raw, KO_abundance, by=c("Bin_name", "KO")) %>%
    distinct()
  
  return(final_table)
}
  



