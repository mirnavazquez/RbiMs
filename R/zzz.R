# zzz.R
.onLoad <- function(...) { }

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    # columnas que aparecen en logs:
    "Abundance", "Abundance_metabolism", "KO", "New", "new2",
    # y otras usadas por dplyr/tidyselect en el paquete:
    "Bin_name", "Percentage", "Pfam", "INTERPRO", "dbCAN_family",
    "rbims_pathway", "rbims_sub_pathway", "Module", "Pathway",
    "Enzyme", "Cycle", "Pathway_cycle", "Detail_cycle", "Genes",
    "Gene_description", "Module_description", "Pathway_description",
    "Scaffold_name", "tmp", "n",
    # dbCAN / MEROPS / helpers (del NOTE de R CMD check)
    "domain_name", "number_of_tools", "Number_of_Tools", "gene_id",
    "dbNamesHMM", "dbNamesdiamond", "dbCAN_names", "dbCAN", "signalp",
    "e_cami", "ecami2", "dbNameseCAMI",
    "after_hyphen", "species_name", "peptidase_code"
  ))
}

