#' @title Read the output of KofamScan/KofamKoala or KAAS
#' @description read_ko() calculates the abundance of each KO within the 
#' bins based on the KofamScan or KofamKoala output.
#' @usage read_ko(data_kofam=NULL, data_kaas=NULL)
#' @param data_kofam a data frame with 5 columns. The contigs are expected to 
#' indicate in their names the bin name and then the scaffold name divided 
#' by an "_", similar to bin_scaffoldXX.
#' @param data_kaas a data frame with 2 columns. The contigs are expected to 
#' indicate in their names the bin name and then the scaffold name divided 
#' by an "_", similar to bin_scaffoldXX.
#' @details This function is part of a package used for the analysis 
#' of bins metabolism.
#' @import dplyr tidyr readr stringr rlang
#' @importFrom utils read.table
#' @examples
#' \dontrun{
#' read_ko("KEGG_bins.txt")
#' }
#' @export
read_ko<-function(data_kofam=NULL, data_kaas=NULL){
  # Kofam_fun -------------------------------------------------------------####
  fun_Kofam<-function(data_kofam){
    # Read data -------------------------------------------------------------####
    KO_raw<-suppressWarnings(
      suppressMessages(
        read_table2(data_kofam, col_names = F) %>%
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
  # Kaas_fun --------------------------------------------------------------####
  fun_kass<-function(data_kaas){
    # Read data -------------------------------------------------------------####
    KAAS<-read.table(data_kaas, sep ="\t", header = F, fill=T)
    
    KAAS_tibble<-as_tibble(KAAS) %>%
      na_if("") %>%
      drop_na() %>%
      rename(Bin_name = .data$V1) %>%
      rename(KO = .data$V2) %>%
      separate(.data$Bin_name, c("Bin_name", "Scaffold_name"), 
               sep = "[_|-][s|S]caffold|-S") %>%
      mutate(Scaffold_name = paste0( "scaffold", .data$Scaffold_name), 
             .data$Scaffold_name)
    # Calculate abundance ---------------------------------------------------####
    KO_abundance<-KAAS_tibble %>%
      group_by(.data$Bin_name) %>%
      distinct() %>%
      count(.data$KO) %>%
      rename(Abundance = .data$n)%>%
      ungroup()
    # Write tibble -----------------------------------------------------------####
    final_table <- left_join(KAAS_tibble, KO_abundance, by=c("Bin_name", "KO")) %>%
      distinct()
    return(final_table)
  }
  # Kofam_Kaas_fun --------------------------------------------------------####
  fun_kofamKaas<-function(data_kofam, data_kaas){
    # Read data KOFAM------------------------------------------------------####
    KO_raw<-suppressWarnings(
      suppressMessages(
        read_table2(data_kofam, col_names = F) %>%
          filter(!str_detect(.data$X1, '#')) %>%
          select(.data$X2, .data$X3) %>%
          rename(Bin_name = .data$X2) %>%
          rename(KO = .data$X3) %>%
          separate(.data$Bin_name, c("Bin_name", "Scaffold_name"),
                   sep = "[_|-][s|S]caffold|-S") %>%
          mutate(Scaffold_name = paste0( "scaffold", .data$Scaffold_name),
                 .data$Scaffold_name)))
    # Read data KAAS ------------------------------------------------------####
    KAAS<-read.table(data_kaas, sep ="\t", header = F, fill=T) %>%
      as_tibble() %>%
      na_if("") %>%
      drop_na() %>%
      rename(Bin_name = .data$V1) %>%
      rename(KO = .data$V2) %>%
      separate(.data$Bin_name, c("Bin_name", "Scaffold_name"), 
               sep = "[_|-][s|S]caffold|-S") %>%
      mutate(Scaffold_name = paste0( "scaffold", .data$Scaffold_name), 
             .data$Scaffold_name)
    # Combine data --------------------------------------------------------####
    combined_data<-bind_rows(KO_raw, KAAS) %>%
      distinct()
    # Calculate abundance -------------------------------------------------####
    KO_abundance<-combined_data %>%
      group_by(.data$Bin_name) %>%
      distinct() %>%
      count(.data$KO) %>%
      rename(Abundance = .data$n)%>%
      ungroup()
    # Write tibble -----------------------------------------------------------####
    final_table <- left_join(combined_data, KO_abundance, by=c("Bin_name", "KO")) %>%
      distinct()
    return(final_table)
  }
  # Evaluating data -------------------------------------------------------####
  if(is.null(data_kofam) == F && is.null(data_kaas) == F){
    final_table<-fun_kofamKaas(data_kofam, data_kaas)
  } else if (is.null(data_kofam) == T && is.null(data_kaas) == F){
    final_table<-fun_kass(data_kaas)
  } else if (is.null(data_kofam) == F && is.null(data_kaas) == T){
    final_table<-fun_Kofam(data_kofam)
  }
  return(final_table)
}