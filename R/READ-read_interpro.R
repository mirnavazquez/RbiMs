#' @title Extract abundance profile of InterProScan output.
#' @description Reads a table object created with InterProScan and generates
#' a profile table of abundance with the hits of the KEGG, PFAM or 
#' INTERPRO databases. The output of KEGG database can be used within 
#' mapping_ko.
#' @usage read_interpro(data_interpro, 
#' database=c("KEGG", "Pfam", "INTERPRO",
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
#' \dontrun{
#' read_interpro(data_interpro="inst/extdata/Interpro_test.tsv", database="INTERPRO", 
#' profile = F)
#' }
#' @export
read_interpro<-function(data_interpro, 
                        database=c("KEGG", "Pfam", "INTERPRO", "TIGRFAM", 
                                   "SUPERFAMILY", "SMART", "SFLD", "ProSiteProfiles",
                                   "ProSitePatterns", "ProDom", "PRINTS", "PIRSF", 
                                   "MobiDBLite","Hamap", "Gene3D", "Coils", "CDD"),
                        profile=TRUE){
  possible_databases<-c("TIGRFAM", "SUPERFAMILY", "SMART", "SFLD", "ProSiteProfiles",
                        "ProSitePatterns", "ProDom", "PRINTS", "PIRSF", "MobiDBLite",
                        "Hamap", "Gene3D", "Coils", "CDD", "Pfam")
  if(database == "KEGG") {
    # Extract KEGG----------------------------------------------####
    table_interpro<-suppressWarnings(
      suppressMessages(read_delim(data_interpro,
                                  delim="\t", col_names = F) %>%
                         drop_na(.data$X12) %>%
                         select(.data$X1, .data$X15) %>%
                         drop_na() %>%
                         distinct() %>%
                         separate_rows(.data$X15, sep="\\|") %>%
                         filter(str_detect(.data$X15, "KEGG")) %>%
                         mutate(X15=str_replace_all(.data$X15, 
                                                    "KEGG: ", "map")) %>%
                         separate(.data$X15, into=c("Pathway", "Enzyme"), 
                                  sep="\\+", 
                                  remove=T, extra ="merge") %>%
                         distinct() %>%
                         rename(Scaffold_full_name = .data$X1)
      ))
    
    # Check enzymes -------------------------------------------------------####
    while(any(str_detect(table_interpro$Enzyme, pattern = "\\+")) == T){
      table_interpro<-table_interpro %>%
        separate_rows(.data$Enzyme, sep="\\+") %>%
        distinct()
    } 
    interpro<- table_interpro %>%
      mutate(Enzyme=str_replace_all(.data$Enzyme, "^", "ec:"))
  } else if (database  %in% possible_databases) {
    # Extract other databases ---------------------------------------------####
    table_interpro_1<-suppressWarnings(
      suppressMessages(read_delim(data_interpro,
                                  delim="\t", col_names = F) %>%
                         drop_na(.data$X12) %>%
                         filter(.data$X4 == database ) %>%
                         select(.data$X1, .data$X5, .data$X6) %>%
                         distinct() %>%
                         rename(Bin_name = .data$X1) %>%
                         rename(!!database := .data$X5) %>%
                         rename(domain_name = .data$X6) 
      ))
    
    table_interpro<-table_interpro_1 %>%
      calc_abundance(analysis = database)
    
    if(isTRUE(profile)){
      interpro<-table_interpro %>%
        select(-.data$Scaffold_name) %>%
        distinct() %>%
        pivot_wider(names_from= "Bin_name", values_from="Abundance", 
                    values_fill = 0)
    } else{
      interpro<-table_interpro
    }
  } else if (database == "INTERPRO") {
    # Extract Interpro ----------------------------------------------------####
    table_interpro_1<-suppressWarnings(
      suppressMessages(read_delim(data_interpro,
                                  delim="\t", col_names = F) %>%
                         drop_na(.data$X12) %>%
                         select(.data$X1, .data$X12, .data$X13) %>%
                         distinct() %>%
                         rename(Bin_name = .data$X1) %>%
                         rename(Interpro = .data$X12) %>%
                         rename(domain_name = .data$X13) 
      ))
    
    table_interpro<-table_interpro_1 %>%
      calc_abundance(analysis = "INTERPRO")
    
    if(isTRUE(profile)){
      interpro<-table_interpro %>%
        select(-.data$Scaffold_name) %>%
        distinct() %>%
        pivot_wider(names_from= "Bin_name", values_from="Abundance", 
                    values_fill = 0)
    } else{
      interpro<-table_interpro
    }
  }
  # Return ----------------------------------------------------------------####
  return(interpro)
}