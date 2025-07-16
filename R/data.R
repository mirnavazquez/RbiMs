#' Main Rbims
#'
#' @format Home made database of intresting pathways
#' @source \url{https://sites.utexas.edu/baker-lab/}
"rbims"

#' Output of read_ko
#' 
#' A dataset containing the output of read_ko.
#'
#' @format A data frame with 6809 rows and four columns.
#' @source \url{https://sites.utexas.edu/baker-lab/}
"ko_bin_table"

#' Output of mapping_ko
#'
#' @description A dataset containing the output of  mapping_ko.
#'
#' @format A data frame containing the pathway profile of six bins.
#' 
#' @source \url{https://sites.utexas.edu/baker-lab/}
"ko_bin_mapp"

#' Metadata
#'
#' @description A dataset containing the metadata of six bins.
#'
#' @format A  data frame of extra data of the bins.
#' 
#' @source \url{https://sites.utexas.edu/baker-lab/}
"metadata"

#' interpro_pfam_profile
#'
#' @description The output of read_interpro for PFAM.
#'
#' @format A  data frame of extra data of the bins.
#' 
#' @source \url{https://sites.utexas.edu/baker-lab/}
"interpro_pfam_profile"

#' interpro_pfam_long
#'
#' @description The output of read_interpro for PFAM.
#'
#' @format A  tiible of extra data of the bins.
#' 
#' @source \url{https://sites.utexas.edu/baker-lab/}
"interpro_pfam_long"

#' interpro_INTERPRO_profile
#'
#' @description The output of read_interpro for Interpro.
#'
#' @format A  tiible of extra data of the bins.
#' 
#' @source \url{https://sites.utexas.edu/baker-lab/}
"interpro_INTERPRO_profile"

#' interpro_INTERPRO_long
#'
#' @description The output of read_interpro for Interpro.
#'
#' @format A  tiible of extra data of the bins.
#' 
#' @source \url{https://sites.utexas.edu/baker-lab/}
"interpro_INTERPRO_long"

#' interpro_KEGG_long
#'
#' @description The output of read_interpro for KEGG
#'
#' @format A  tiible of extra data of the bins.
#' 
#' @source \url{https://sites.utexas.edu/baker-lab/}
"interpro_KEGG_long"

#' interpro_map
#'
#' @description The output of mapping for KEGG from read_interpro.
#'
#' @format A tiible of extra data of the bins.
#' 
#' @source \url{https://sites.utexas.edu/baker-lab/}
"interpro_map"
#' KO abundances predicted with PICRUSt2
#'
#' This dataset contains KO (KEGG Orthology) predicted abundances from PICRUSt2 output.
#' @format A data frame with X rows and Y columns (describe briefly).
#' @source Generated with PICRUSt2 using X sample set.
#' @docType data
#' @name KO_picrust2
"KO_picrust2"

#' EC abundances predicted with PICRUSt2
#'
#' Dataset with EC number abundances predicted from functional profiling.
#' @format A data frame with X rows and Y columns.
#' @source Generated from PICRUSt2.
#' @docType data
#' @name EC_picrust2
"EC_picrust2"

#' Pathway-level predictions with PICRUSt2
#'
#' Pathway predictions from MetaCyc using PICRUSt2 output.
#' @format A data frame with X rows and Y columns.
#' @source PICRUSt2 results.
#' @docType data
#' @name pathway_picrust2
"pathway_picrust2"

