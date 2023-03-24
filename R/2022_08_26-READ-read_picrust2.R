#' @title Read the output of PICRUSt2 
#' @description read a table created in PICRUSt2 and gives format to downstream analysis
#' @usage read_picrust2(data_picrust2, profile=TRUE, write=FALSE, 
#' database= c("KO", "EC", "pathway"))
#' @param data_picrust2 a table, output of PICRUSt2 on tsv format.
#' @param database a character indicating for which database do you want to
#' get the abundance profile. Valid options are "KO", "EC" or "pathway".
#' @param profile a logical value indicating if you want to print a profile 
#' or not.
#' @param write  a logical value indicating to save the data imported 
#' as a formatted table with .tsv extension with a time stamp
#' @details This function is part of a package used for the analysis of bins 
#' metabolism.
#' @import tibble dplyr stringr tidyr janitor rlang
#' @examples
#' \dontrun{
#' read_picrust2(data_picrust2="inst/extdata/pred_metagenome_unstrat.tsv", 
#' database="KO", profile = F, write=F)
#' }
#' @export
read_picrust2<-function(data_picrust2, 
                        profile=TRUE, write=FALSE,
                        database= c("KO", "EC", "pathway")){
# Extract functions or gene----------------------------------------------####
    table_picrust<-suppressWarnings(
      suppressMessages(read_delim(data_picrust2,
                                  delim="\t", 
                                  col_names = T)  )) 
# Choosing database----------------------------------------------####
    
if(database=="KO"){
          colnames(table_picrust)[1] <- "KO"}
if(database=="EC"){
      colnames(table_picrust)[1] <- "EC"}   
if(database=="pathway"){
     colnames(table_picrust)[1] <- "pathway"}   
    
# Profile or not ----------------------------------------------####
    
 if(isFALSE(profile)){
      picrust<-table_picrust%>% 
        pivot_longer(-1,names_to= "Bin_name", 
                     values_to="Abundance")
      
    } else{
      picrust<-table_picrust
    }

# Write data or not --------------------------------------------------------------####

if(isTRUE(write)){
  write_tsv(picrust, paste0("picrust_output_", format(Sys.time(), "%b_%d_%X"), ".tsv"))
}
else{
  return(picrust)
}

# Return ----------------------------------------------------------------####
return(picrust)
}


