#' @title Extract abundance profile of InterProScan output.
#' @description Reads a table object created with InterProScan and generates
#' a profile table of abundance with the hits of the KEGG, Pfam or 
#' INTERPRO databases. The output of KEGG database can be used within 
#' mapping_ko.
#' @usage read_interpro (data_interpro, 
#' database = c("KEGG", "Pfam", "INTERPRO",
#' "TIGRFAM", "SUPERFAMILY", "SMART", "SFLD", "ProSiteProfiles",
#' "ProSitePatterns", "ProDom", "PRINTS", "PIRSF", 
#' "MobiDBLite","Hamap", "Gene3D", "Coils", "CDD"), profile = TRUE)
#' @param data_interpro a table, output of InterProScan on tsv format.
#' InterProScan should have been run with -pa option to be able to use the 
#' KEGG option, in the database argument.
#' @param database a character indicating for which database do you want to
#' get the abundance profile. Valid options are "KEGG", "PFAM" or "INTERPRO".
#' @param profile a logical value indicating if you want to print a profile 
#' or not. This option is valid for "PFAM" and "INTERPRO" database. 
#' @details This function is part of a package used for the analysis of bins 
#' metabolism.
#' @import tibble dplyr stringr tidyr janitor rlang
#' @examples
#' read_interpro(data_interpro = "inst/extdata/Interpro_test.tsv", 
#'               database = "INTERPRO", profile = FALSE)

#' @export
read_interpro<-function(data_interpro = NULL, 
                        database=c("KEGG", "Pfam", "INTERPRO", "TIGRFAM", 
                                   "SUPERFAMILY", "SMART", "SFLD", "ProSiteProfiles",
                                   "ProSitePatterns", "ProDom", "PRINTS", "PIRSF", 
                                   "MobiDBLite","Hamap", "Gene3D", "Coils", "CDD"),
                        profile=TRUE){
  ruta_interpro <- data_interpro
  possible_databases<-c("TIGRFAM", "SUPERFAMILY", "SMART", "SFLD", "ProSiteProfiles",
                        "ProSitePatterns", "ProDom", "PRINTS", "PIRSF", "MobiDBLite",
                        "Hamap", "Gene3D", "Coils", "CDD", "Pfam")
  if(database == "PFAM"){
    stop("P letter must be in capital followed by lower-case letters e.g 'Pfam'")
  }
  if(database == "KEGG") {
    # Extract KEGG----------------------------------------------#### Aquí también 
    lapply_read_delim_bind_rows <- function(path, pattern = "*.tsv"){
      if (file.exists(path) && grepl("\\.tsv$", path)) {
        return(read_tsv(path, col_names = FALSE))
      }
      files <- list.files(path, pattern = "*.tsv", full.names = TRUE, recursive = TRUE)
      lapply(files, read_tsv, col_names = FALSE) %>%
        bind_rows()
    }
    
    final_files <- suppressWarnings(lapply_read_delim_bind_rows(ruta_interpro)) %>%
        suppressWarnings() %>%
        suppressMessages() %>%
        drop_na(.data$X12) %>%
        select(.data$X1, .data$X15) %>%
        drop_na() %>%
        distinct() %>%
        separate_rows(.data$X15, sep="\\|") %>%
        filter(str_detect(.data$X15, "KEGG")) %>%
        mutate(X15=str_replace_all(.data$X15, "KEGG: ", "map")) %>%
        separate(.data$X15, into=c("Pathway", "Enzyme"), sep="\\+", 
                                  remove=T, extra ="merge") %>%
        distinct() %>%
        rename(Scaffold_full_name = .data$X1)
    table_interpro <- final_files 

    
    # Check enzymes -------------------------------------------------------####
    while(any(str_detect(table_interpro$Enzyme, pattern = "\\+")) == T){
      table_interpro<-table_interpro %>%
        separate_rows(.data$Enzyme, sep="\\+") %>%
        distinct()
    } 
    interpro<- table_interpro %>%
      mutate(Enzyme=str_replace_all(.data$Enzyme, "^", "ec:"))
  } else if (database  %in% possible_databases) { 
    
    # Extract other databases ---------------------------------------------#### Mas o menos por aqui
    lapply_read_delim_bind_rows <- function(path, pattern = "*.tsv"){
      if (file.exists(path) && grepl("\\.tsv$", path)) {
        return(read_tsv(path, col_names = FALSE))
      }
      files <- list.files(path, pattern = "*.tsv", full.names = TRUE, recursive = TRUE)
      lapply(files, read_tsv, col_names = FALSE) %>%
        bind_rows()
    }
    
    final_files <- suppressWarnings(lapply_read_delim_bind_rows(ruta_interpro)) %>%
        suppressWarnings() %>%
        suppressMessages() %>%
        drop_na(.data$X12) %>%
        filter(.data$X4 == database ) %>%
        select(.data$X1, .data$X5, .data$X6) %>%
        distinct() %>%
        rename(Bin_name = .data$X1) %>%
        rename(!!database := .data$X5) %>%
        rename(domain_name = .data$X6) %>%
        calc_abundance(analysis = database, col_rename = database)
    table_interpro <- final_files 

    if(isTRUE(profile)){
      interpro<-table_interpro %>%
        select(-.data$Scaffold_name) %>%
        distinct() %>%
        pivot_wider(names_from= "Bin_name", values_from="Abundance", 
                    values_fill = 0)
    } else{
      
      interpro <- table_interpro
    }
  } else if (database == "INTERPRO") {
    # Extract Interpro ----------------------------------------------------####
    lapply_read_delim_bind_rows <- function(path, pattern = "*.tsv"){
      if (file.exists(path) && grepl("\\.tsv$", path)) {
        return(read_tsv(path, col_names = FALSE))
      }
      files <- list.files(path, pattern = "*.tsv", full.names = TRUE, recursive = TRUE)
      lapply(files, read_tsv, col_names = FALSE) %>%
        bind_rows()
    }
    final_files <- suppressWarnings(lapply_read_delim_bind_rows(ruta_interpro)) %>%
        suppressWarnings() %>%
        suppressMessages() %>%
        drop_na(.data$X12) %>%
        select(.data$X1, .data$X12, .data$X13) %>%
        distinct() %>%
        rename(Bin_name = .data$X1) %>%
        rename(Interpro = .data$X12) %>%
        rename(domain_name = .data$X13) %>%
        calc_abundance(analysis = "INTERPRO", col_rename = "INTERPRO")
     table_interpro_1 <-final_files

    
    if(isTRUE(profile)){
      interpro<-table_interpro_1 %>%
        select(-.data$Scaffold_name) %>%
        distinct() %>%
        pivot_wider(names_from= "Bin_name", values_from = "Abundance", 
                    values_fill = 0)
    } else{
      
      interpro <- table_interpro_1
    }
  }
  # Return ----------------------------------------------------------------####
  return(interpro)
}
