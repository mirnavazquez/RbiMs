#' @title Read the output of dbCAN and extract abundance profile 
#' @description read_dbcan3 calculates the abundance of each Gene within the 
#' bins based on the dbCAN output from run_dbcan3 V3.0.6. 
#' @usage read_dbcan3(dbcan_path, write=FALSE, profile=TRUE)
#' @param dbcan_path a path where dbCAN output data are. They 
#' should have the extension overview.txt and all files in the path are the ones that
#' need to be read. Output data should have 6 columns with the bin names 
#' followed by the Genes obtained in every algorithm (HMMER,Hotpep,DIAMOND), 
#' column 'Signalp' indcating if a Peptide signal is found and a column 
#' '#ofTools" indicating the number of algorithms that found this Gene. 
#' @param profile a logical value indicating if you want to print a profile 
#' or not.
#' @param write  a logical value indicating to save the data imported 
#' as a formatted table with .tsv extension with a time stamp and it will be 
#' located in your current workin directory
#' @import dplyr tidyr readr stringr rlang tidyselect utils 
#' @examples
#' \dontrun{
#' read_dbcan3("C:/Users/bins/")
#' }
#' @export

read_dbcan3<-function(dbcan_path, write=FALSE, profile=TRUE){
  ruta_dbcan<-dbcan_path
  
  # Load all the data tables results ---------------------------------------####
  lapply_read_delim_bind_rows <- function(path, pattern = "overview.txt"){
    files = list.files(path, pattern ="overview.txt", full.names = TRUE, recursive = T)
    lapply(files, read.delim, check.names=F) %>% 
      bind_rows()
  }
  dbcan_df<-suppressWarnings(lapply_read_delim_bind_rows(ruta_dbcan)) 
  
  
  # Reading data ----------------------------------------------------------####
  dbcan_df_format<- suppressWarnings(
    dbcan_df %>%
      clean_names() %>%
        rename(Number_of_Tools = number_of_tools)  %>%
        filter(Number_of_Tools > 1)  %>%
        rename(Bin_name = gene_id) %>%
        mutate( hmmer2=str_replace_all(.data$hmmer, "[[:punct:]]",  "\t")) %>%
        separate(.data$hmmer2, c("dbNamesHMM"),  sep="\t") %>%
        mutate(diamond2=str_replace_all(.data$diamond, "[[:punct:]]", "\t")) %>%
        separate(.data$diamond2, c("dbNamesdiamond"), sep="\t"))

  ###################################### e_cami ##################################  
  if (is.null(dbcan_df_format$e_cami)) { 
    df_format <- dbcan_df_format %>%
      unite("dbCAN_names", dbNamesHMM, dbNamesdiamond, sep = "_", remove = FALSE) %>%  
      mutate(dbCAN_names = str_replace_all(dbCAN_names, "^_", "")) %>%
      separate(dbCAN_names, c("dbCAN_names"), sep = "_") %>%
      mutate(dbCAN = case_when(
        str_detect(dbCAN_names, "CBM") ~ "carbohydrate-binding module [CBM]",
        str_detect(dbCAN_names, "CE") ~ "carbohydrate esterases [CEs]",
        str_detect(dbCAN_names, "GH") ~ "glycoside hydrolases [GHs]",
        str_detect(dbCAN_names, "GT") ~ "glycosyltransferases [GTs]",
        str_detect(dbCAN_names, "PL") ~ "polysaccharide lyases [PLs]",
        str_detect(dbCAN_names, "AA") ~ "auxiliary activities [AAs]"
      )) %>%
      mutate_if(is.character, str_trim) %>%
      select(Bin_name, dbCAN_names, domain_name = dbCAN, signalp) %>%
      calc_abundance(analysis = "dbCAN", col_rename = "dbCAN_names") %>%
      select(-Scaffold_name, dbCAN_names) %>%
      group_by(Bin_name, dbCAN_family = dbCAN_names, domain_name, signalp) %>%
      summarise_if(is.numeric, sum)
  } else {
    df_format <- dbcan_df_format %>%
      mutate(ecami2 = str_replace_all(e_cami, "[[:punct:]]", "\t")) %>%
      separate(ecami2, c("dbNameseCAMI"), sep = "\t") %>%
      unite("dbCAN_names", dbNamesHMM, dbNamesdiamond, dbNameseCAMI, sep = "_", remove = FALSE) %>%
      mutate(dbCAN_names = str_replace_all(dbCAN_names, "^_", "")) %>%
      separate(dbCAN_names, c("dbCAN_names"), sep = "_") %>%
      mutate(dbCAN = case_when(
        str_detect(dbCAN_names, "CBM") ~ "carbohydrate-binding module [CBM]",
        str_detect(dbCAN_names, "CE") ~ "carbohydrate esterases [CEs]",
        str_detect(dbCAN_names, "GH") ~ "glycoside hydrolases [GHs]",
        str_detect(dbCAN_names, "GT") ~ "glycosyltransferases [GTs]",
        str_detect(dbCAN_names, "PL") ~ "polysaccharide lyases [PLs]",
        str_detect(dbCAN_names, "AA") ~ "auxiliary activities [AAs]"
      )) %>%
      mutate_if(is.character, str_trim) %>%
      select(Bin_name, dbCAN_names, domain_name = dbCAN, signalp) %>%
      calc_abundance(analysis = "dbCAN", col_rename = "dbCAN_names") %>%
      select(-Scaffold_name, dbCAN_names) %>%
      group_by(Bin_name, dbCAN_family = dbCAN_names, domain_name, signalp) %>%
      summarise_if(is.numeric, sum)
  }

  # Reformating data ----------------------------------------------------------####
  
  dbcan_df_reformat <-df_format %>%
    dplyr::select(-.data$signalp, .data$dbCAN_family) %>%
    group_by(.data$Bin_name, .data$dbCAN_family, .data$domain_name) %>% 
    summarize_if(is.numeric, sum) %>%   
    pivot_wider(names_from = "Bin_name", values_from = "Abundance") %>% 
    ungroup() %>% 
    mutate_if(is.numeric, ~replace(., is.na(.), 0)) 
  
  # Messages --------------------------------------------------------------####
  initial<-dim(dbcan_df)
  final<-dim(dbcan_df_format %>%
               filter(Number_of_Tools >1))
  signals<- df_format %>% 
    group_by(.data$signalp) %>% 
    count() 
  signals2<- df_format %>% 
    group_by(.data$signalp) %>% 
    count() 
  
  
  print(paste0("Input Genes = " , initial[1]))
  print(paste0("Remained Genes after filtering = " , final[1]))
  print(paste0("Percentage of genes remained = " , round(final[1]/initial[1]*100), "%"))
  print(paste0("Number of genes with signals = " , sum(signals[-1,]$n)))
  print(paste0("Number of genes with signals that passed filtering = " , sum(signals2[-1,]$n)))
  
  # Profile or not --------------------------------------------------------------####
  
  if(isTRUE(profile)){
    output<-dbcan_df_reformat
  } else{
    output<-df_format 
  }
  
  
  # Write data or not --------------------------------------------------------------####
  
  if(isTRUE(write)){
    write_tsv(output, paste0("dbcan_output_", format(Sys.time(), "%b_%d_%X"), ".tsv"))
  }
  else{
    return(output)
  }
  
  # Return ----------------------------------------------------------------####
  return(output)
  
}
