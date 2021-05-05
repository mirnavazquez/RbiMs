#' @title Read the output of KofamScan/KofamKoala or KAAS.
#' @description read_ko calculates the abundance of each KO within the 
#' bins based on the KofamScan or KofamKoala output.
#' @usage read_ko(data_kofam=NULL, data_kaas=NULL, data_interpro=NULL)
#' @param data_kofam a data frame with 5 columns. Contigs are expected to 
#' indicate in their names the bin name follow by the scaffold name 
#' divided by an "_": bin_scaffoldXX.
#' @param data_kaas a data frame with 2 columns. Contigs are expected to 
#' indicate in their names the bin name follow by the scaffold name 
#' divided by an "_": bin_scaffoldXX. 
#' @param data_interpro a data frame output of read_interpro. This
#' argument is used within mapping_KO.
#' @details This function is part of a package used for the analysis 
#' of bins metabolism.
#' @import dplyr tidyr readr stringr rlang
#' @importFrom utils read.table
#' @examples
#' \dontrun{
#' read_ko("KEGG_bins.txt")
#' }
#' @export
read_ko<-function(data_kofam=NULL, 
                  data_kaas=NULL, 
                  data_interpro=NULL){
  # Kofam_fun -------------------------------------------------------------####
  if( is.null(data_kofam) == F && is.null(data_kaas) == F || 
      is.null(data_kofam) == F){
    table_Kofam<-suppressWarnings(
      suppressMessages(
        read_table2(data_kofam, col_names = F) %>%
          filter(!str_detect(.data$X1, '#')) %>%
          select(.data$X2, .data$X3) %>%
          rename(Bin_name = .data$X2) %>%
          rename(KO = .data$X3)
      ))
  }
  
  # Kaas_fun --------------------------------------------------------------####
  if(is.null(data_kofam) == F && is.null(data_kaas) == F || 
     is.null(data_kaas) == F){
    table_KAAS<-read.table(data_kaas, sep ="\t", header = F, fill=T) %>%
      as_tibble() %>%
      na_if("") %>%
      drop_na() %>%
      rename(Bin_name = .data$V1) %>%
      rename(KO = .data$V2)
  }
  # Interpro --------------------------------------------------------------####
  if (is.null(data_interpro) == F){
    tabla_to_print <- data_interpro %>%
      select(.data$Scaffold_full_name, .data$KO) %>%
      drop_na() %>%
      rename(Bin_name = .data$Scaffold_full_name) %>%
      calc_abundance(analysis="KEGG")
  }
  # Calc ------------------------------------------------------------------####
  if(is.null(data_kofam) == F && is.null(data_kaas) == F ){
    tabla_to_print<-bind_rows(table_Kofam, table_KAAS) %>%
      distinct() %>%
      calc_abundance(analysis="KEGG")
  }
  if( is.null(data_kofam) == F){
    tabla_to_print<-calc_abundance(table_Kofam, analysis="KEGG")
  }
  if( is.null(data_kaas) == F){
    tabla_to_print<-calc_abundance(table_KAAS, analysis="KEGG")
  }
  return(tabla_to_print)
}
