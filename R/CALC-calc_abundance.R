#' @title Calculate Abundance
#' @description This function calculates the abundance of a certain feature, 
#' based on the number of times it appears in a dataset.
#' @usage calc_abundance(dataset, analysis=c("KEGG", "PFAM", "INTERPRO", "dbCAN"))
#' @param dataset A data frame with two columns: the gene name
#'  and the associated ID.
#' @param analysis A character vector indicating which type of analysis 
#' is being performed.
#' @param col_rename  A character vector indicating which type of analysis 
#' is being performed, same as analysis.
#' @details This function calculates the abundance of a feature by summing 
#' the number of times it appears in the dataset. The formula used is Σx, 
#' where x represents the number of times the feature appears.
#' @import tibble dplyr stringr tidyr rlang
#' @return A data frame with four columns: the gene ID, the genome name, 
#' the KEGG ID (or other feature ID), and the number of times the feature appears 
#' in the genome.
#' @noRd
calc_abundance <- function(dataset, 
                           analysis = c("KEGG", "Pfam", "INTERPRO", "TIGRFAM", 
                                        "SUPERFAMILY", "SMART", "SFLD",
                                        "ProSiteProfiles", "ProSitePatterns", 
                                        "ProDom", "PRINTS", "PIRSF", 
                                        "MobiDBLite", "Hamap", "Gene3D", 
                                        "Coils", "CDD", "dbCAN"), 
                           col_rename = NULL) {
  # Asign colum names -----------------------------------------------------####
  col_analysis <- c(KEGG = "KO", Pfam = "Pfam", INTERPRO = "Interpro", 
                    dbCAN = "dbCAN_names", TIGRFAM = "TIGRFAM",
                    SUPERFAMILY = "SUPERFAMILY", SMART = "SMART", SFLD = "SFLD",
                    ProSiteProfiles = "ProSiteProfiles", 
                    ProSitePatterns = "ProSitePatterns", ProDom = "ProDom", 
                    PRINTS = "PRINTS", PIRSF = "PIRSF", 
                    MobiDBLite = "MobiDBLite",  Hamap = "Hamap", 
                    Gene3D = "Gene3D",  Coils = "Coils", CDD = "CDD")
  if (!analysis %in% names(col_analysis)) stop("Unknown analysis")
  # Read table ------------------------------------------------------------####
  KO_raw <- dataset %>%
    separate(Bin_name, c("Bin_name", "Scaffold_name"), 
             sep = "[_|-][s|S]caffold") %>%
    mutate(Scaffold_name = paste0("scaffold", Scaffold_name),
           Scaffold_name) %>%
    unite("Scaffold_name", c("Bin_name", "Scaffold_name"), remove = FALSE)
  # Selecting analysis ----------------------------------------------------####
  KO_raw <- KO_raw %>%
    rename(tmp = .data[[col_analysis[analysis]]]) %>%
    distinct()
  # Calculate abundance ---------------------------------------------------####
  KO_abundance <- KO_raw %>%
    group_by(Bin_name) %>%
    distinct() %>%
    count(tmp) %>%
    rename(Abundance = n) %>%
    ungroup()
  # Write tibble -----------------------------------------------------------####
  final_table_1 <- left_join(KO_raw, KO_abundance, 
                             by = c("Bin_name", "tmp")) %>%
    distinct()
  
  if (!is.null(col_rename)) {
    final_table <- final_table_1 %>% rename_with(~ col_rename, tmp)
  } else {
    final_table <- final_table_1
  }
  
  return(final_table)
}

