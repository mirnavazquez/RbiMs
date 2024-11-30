#' @title Read the output of merops and extract abundance profile 
#' @description read_merops calculates the abundance of each protease within the 
#' bins based on the merops output from run_merops. 
#' @usage read_merops(merops_path, write=FALSE, profile=TRUE)
#' @param merops_path a path where Merops output data are. They 
#' should have the extension .txt and all files in the path are the ones that
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

read_merops <- function(merops_path, write = FALSE, profile = TRUE) {
  ruta_merops <- merops_path
  
  # Load all the data tables results ---------------------------------------####
  lapply_read_delim_bind_rows <- function(path, pattern = ".txt") {
    files <- list.files(path, pattern = ".txt", full.names = TRUE, recursive = TRUE)
    lapply(files, read.delim, header = F, check.names = F) %>% 
      bind_rows()
  }
  merops_df <- suppressWarnings(lapply_read_delim_bind_rows(ruta_merops)) 
  

  # Reading data ----------------------------------------------------------####
  merops_df_format <- suppressWarnings(
    merops_df %>%
      filter(.data$V4 > 90) %>%
      rename(Bin_name = .data$V1)) %>%
      mutate(
        after_hyphen = str_extract(.data$V3, "(?<=-\\s).*?(?= \\()"),
        species_name = str_extract(.data$V3, "(?<=\\().*?(?=\\))"),
        peptidase_code = str_extract(.data$V3, "(?<=\\[).*?(?=\\])"),
      MEROPS_names = paste(after_hyphen, species_name, peptidase_code, sep = " | ")
      ) %>%
      rename(domain_name = .data$V2)%>%
      calc_abundance(analysis = "MEROPS", col_rename = "MEROPS_names") %>%
      dplyr::select(-.data$Scaffold_name, .data$MEROPS_names) %>% 
      group_by(.data$Bin_name, MEROPS_family = .data$MEROPS_names, 
             .data$domain_name) %>% 
      summarise_if(is.numeric, sum)
  
  # Reformating data ----------------------------------------------------------####
  
  merops_df_reformat <-merops_df_format %>%
    select(-.data$V4, -.data$V5, -.data$V6) %>%
    group_by(.data$Bin_name, .data$MEROPS_family, .data$domain_name) %>% 
    summarize_if(is.numeric, sum) %>%   
    pivot_wider(names_from = "Bin_name", values_from = "Abundance") %>% 
    ungroup() %>% 
    mutate_if(is.numeric, ~replace(., is.na(.), 0)) 
  
  # Profile or not --------------------------------------------------------------####
  
  if(isTRUE(profile)){
    output<-merops_df_reformat
  } else{
    output<-merops_df_format%>%
      select(-.data$V4, -.data$V5, -.data$V6)
  }
  
  # Write data or not --------------------------------------------------------------####
  
  if(isTRUE(write)){
    write_tsv(output, paste0("merops_output_", format(Sys.time(), "%b_%d_%X"), ".tsv"))
  }
  else{
    return(output)
  }
  
  # Return ----------------------------------------------------------------####
  return(output)
  
}
  

  