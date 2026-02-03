#' @title Metabolism PCA
#' @description Identifies the most important KEGG pathways or protein 
#' domains in the entire database and returns a profile of those with 
#' the highest contributions to PCA dimensions.
#' @usage get_subset_pca(tibble_rbims, cos2_val = NULL, 
#' analysis = c("KEGG", "Pfam", "INTERPRO", "dbCAN", "MEROPS"))
#' @param tibble_rbims A tibble created with functions such as 
#' read_interpro(), mapping_ko(), or get_subset_*().
#' @param cos2_val Numeric value between 0 and 1 indicating the cutoff 
#' for contribution. Defaults to 0.98. See \link[factoextra]{get_pca}.
#' @param analysis Character string specifying the annotation database. 
#' Options: "KEGG", "Pfam", "INTERPRO", "dbCAN", "MEROPS".
#' @details This function is part of the rbims package for analyzing 
#' metabolic potential in metagenome-assembled genomes (MAGs).
#' @import dplyr factoextra rlang tibble
#' @examples
#' # get_subset_pca(ko_bin_mapp, analysis = "KEGG")
#' @export
get_subset_pca <- function(tibble_rbims,
                           cos2_val = NULL,
                           analysis = c("KEGG", "Pfam", "INTERPRO", "dbCAN", "MEROPS")) {
  
  # Error handling -------------------------------------------------------------
  if (analysis == "PFAM") {
    stop("Use 'Pfam' with a capital P followed by lowercase letters.")
  }
  
  # Select data ----------------------------------------------------------------
  if (analysis == "KEGG") {
    data_to_select <- c("Module", "Module_description", "Pathway", 
                        "Pathway_description", "Genes", "Gene_description", 
                        "Enzyme", "Cycle", "Pathway_cycle", "Detail_cycle", 
                        "rbims_pathway", "rbims_sub_pathway")
  } else if (analysis %in% c("Pfam", "INTERPRO", "dbCAN", "MEROPS")) {
    data_to_select <- "domain_name"
  }
  
  # Rename columns -------------------------------------------------------------
  if (analysis == "Pfam") {
    tibble_rbims <- tibble_rbims %>% 
      rename(tmp = .data$Pfam)
  } else if (analysis == "KEGG") {
    tibble_rbims <- tibble_rbims %>% 
      rename(tmp = .data$KO)
  } else if (analysis == "INTERPRO") {
    tibble_rbims <- tibble_rbims %>% 
      rename(tmp = .data$INTERPRO)
  } else if (analysis == "dbCAN") {
    tibble_rbims <- tibble_rbims %>% 
      rename(tmp = .data$dbCAN_family)
  } else if (analysis == "MEROPS") {
    tibble_rbims <- tibble_rbims %>% 
      rename(tmp = .data$MEROPS_family)
  }
  
  # Ensure unique rownames -----------------------------------------------------
  tibble_rbims$tmp <- make.unique(as.character(tibble_rbims$tmp))
  
  wide_ko <- tibble_rbims %>%
    dplyr::select(-all_of(data_to_select)) %>%
    dplyr::distinct() %>%
    tibble::column_to_rownames("tmp")
  
  # Run PCA --------------------------------------------------------------------
  if (!requireNamespace("FactoMineR", quietly = TRUE)) {
    install.packages("FactoMineR")
  }
  df_pca <- FactoMineR::PCA(wide_ko, scale.unit = TRUE, graph = FALSE)
  contribution <- as.data.frame(df_pca$ind$cos2)
  
  # Warnings and messages ------------------------------------------------------
  if (all(contribution$Dim.1 <= 0.97)) {
    warning("Contribution of the first dimension is <= 0.97. ",
            "If the output is empty, try lowering cos2_val.")
  } else {
    message("Contribution of the first dimension is > 0.97. ",
            "If the output is empty, try adjusting cos2_val.")
  }
  
  
  # Default cutoff -------------------------------------------------------------
  if (is.null(cos2_val)) {
    cos2_val <- 0.98
  }
  
  # Extract features with high contribution ------------------------------------
  subset1 <- rownames(subset(contribution, contribution$Dim.1 >= cos2_val))
  subset2 <- rownames(subset(contribution, contribution$Dim.2 >= cos2_val))
  selected <- unique(c(subset1, subset2))
  
  # Filter tibble and rename back ----------------------------------------------
  final_table <- tibble_rbims %>% filter(.data$tmp %in% selected)
  
  if (analysis == "Pfam") {
    final_table <- final_table %>% rename(Pfam = .data$tmp)
  } else if (analysis == "KEGG") {
    final_table <- final_table %>% rename(KO = .data$tmp)
  } else if (analysis == "INTERPRO") {
    final_table <- final_table %>% rename(INTERPRO = .data$tmp)
  } else if (analysis == "dbCAN") {
    final_table <- final_table %>% rename(dbCAN_family = .data$tmp)
  } else if (analysis == "MEROPS") {
    final_table <- final_table %>% rename(MEROPS_family = .data$tmp)
  }
  
  return(final_table)
}
