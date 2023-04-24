#' @title Read the output of KofamScan/KofamKoala or KAAS.
#' @description read_ko_shiny calculates the abundance of each KO within the 
#' bins based on the KofamScan or KofamKoala output.
#' @usage read_ko_shiny(data_kofam=NULL, data_kaas=NULL, data_interpro=NULL)
#' @param data_kofam a file output of KofamScan/KofamKoala. They 
#' should have the extension .txt and all files in the path are the ones that
#' need to be read. Output data should have 5 columns with the bin names 
#' followed by the scaffold name divided by a '-' or '_': bin_scaffoldXX.
#' @param data_kaas a data frame with 2 columns. Contigs are expected to 
#' indicate in their names the bin name followed by the scaffold name 
#' divided by a '-' or '_': bin_scaffoldXX. 
#' @param data_interpro a data frame output of read_interpro. This
#' argument is used within mapping_KO.
#' @details This function is part of a package used for the analysis 
#' of bins metabolism.
#' @import dplyr tidyr readr stringr rlang
#' @importFrom utils read.table
#' @importFrom purrr map_dfr 
#' @examples
#' \dontrun{
#' read_ko_shiny("Kofam_output")
#' }
#' @export
read_ko_shiny<-function(data_kofam=NULL, 
                  data_kaas=NULL, 
                  data_interpro=NULL){
  # Kofam_fun -------------------------------------------------------------####
  if( is.null(data_kofam) == F && is.null(data_kaas) == F || 
      is.null(data_kofam) == F){
    table_Kofam<-suppressWarnings(
      suppressMessages( 
        read_table2(data_kofam, col_names = F) %>%
          filter(str_detect(.data$X1, '\\*')) %>%
          select(.data$X2, .data$X3) %>%
          rename(Bin_name = .data$X2) %>%
          rename(KO = .data$X3)
      ))
    if (length(grep("[S|s]caffold", table_Kofam$Bin_name, invert = TRUE)) != 0){
      stop("Must label scaffolds with the name 'Scaffold' or 'scaffold' after 
  the bin name followed by a '-' or '_'.")
    }
    before<- sub("scaffold.*", "",table_Kofam$Bin_name)
    extract<- substr(before, nchar(before)-1+1, nchar(before))
    if (all(str_detect(extract, "[-|_]"))){
    } else{
      stop("Bin name and scaffold is not separated by '-' or '_'.")
    }
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
    if (length(grep("[S|s]caffold", table_Kofam$Bin_name, invert = TRUE)) != 0){
      stop("Must label scaffolds with the name 'Scaffold' or 'scaffold' after 
  the bin name followed by a '-' or '_'.")
    }
    before<- sub("scaffold.*", "",table_Kofam$Bin_name)
    extract<- substr(before, nchar(before)-1+1, nchar(before))
    if (all(str_detect(extract, "[-|_]"))){
    } else{
      stop("Bin name and scaffold is not separated by '-' or '_'.")
    }
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
