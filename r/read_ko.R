read_ko<-function(ko_result){
  read_table2(ko_result, col_names = F) %>%
    filter(!str_detect(X1, '#')) %>%
    select(X2, X3) %>%
    rename(Bin_name = X2) %>%
    rename(ko = X3) %>%
    separate(Bin_name, c("Bin_name", "Scaffold_name"), sep = "[_|-][s|S]caffold") %>%
    mutate(Scaffold_name = paste0( "scaffold", Scaffold_name), Scaffold_name)
  }
  
read_ko("data/pruena.txt")


