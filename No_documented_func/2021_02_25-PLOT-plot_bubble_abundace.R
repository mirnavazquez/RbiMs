library(ggplot2)
library(dplyr)

plot_bubble_abundance<-function(tabble_ko,
                                other_data,
                                x_axis,
                                y_axis,
                                size_bubble, 
                                metadata_feature){
  ############################ quoting ##############################
  met_features <- enquo(y_axis)
  met_features_x <- as_label(met_features)
  other_features_x <- as_label(metadata_feature)
  ######################## Parsing the table ########################
    kegg_longer<-tabble_ko %>%
    pivot_longer(cols = -c(Module, Module_description, Pathway, 
                           Pathway_description, Genes, 
                           Gene_description, Enzyme, KO),  
                 names_to = "Bin_name", 
                 values_to = "Abundance") %>%
    distinct() %>%
    mutate_at('Abundance', as.integer) %>%
    select(all_of(met_features_x), Bin_name, Abundance) %>%
    distinct() %>%
    drop_na() %>%
    left_join(other_data, by="Bin_name") %>%
    mutate(Abundance = case_when(
      Abundance == 0 ~ NA_integer_,
      TRUE ~ as.integer(Abundance)
    ))
  ############################## Plot #############################
  ggplot2::ggplot(kegg_longer,
                  aes_string(x=x_axis, 
                             y=met_features, 
                            size=size_bubble,
                            color=metadata_feature)) +
    geom_point(alpha=0.5) +
    scale_size(range = c(1,10))+
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45,
                                     size=8,
                                     hjust = 1,
                                     vjust = 1),
          axis.text.y = element_text(size=7))
}

plot_bubble_abundance(Kegg_subset, 
                      metadata, 
                      "Bin_name",
                       Module,
                      "Abundance",
                      "Clades")
