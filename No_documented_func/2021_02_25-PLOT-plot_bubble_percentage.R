library(ggplot2)
library(dplyr)

plot_bubble_percentage<-function(tabble_ko,
                                other_data, 
                                x_axis,
                                y_axis,
                                size_bubble, 
                                metadata_feature){
  ############################ quoting ##############################
  y_axis_features <- enquo(y_axis)
  met_features_y <- as_label(y_axis_features)
  other_features_x <- as_label(metadata_feature)
  #################### Transform from wide to long ####################
  Kegg_long<- tabble_ko %>%
  pivot_longer(cols = -c(Module, Module_description, Pathway, 
            Pathway_description, Genes, 
            Gene_description, Enzyme, KO), values_to = "Abundance",
            names_to="Bin_name") %>%
  distinct()
  ######### Count the total number of genes per metabolism ########
  metabolism_counts_totals<-Kegg_long %>%
    select(all_of(met_features_y), KO) %>%
    distinct() %>%
    count(!!y_axis_features, sort=T) %>%
    rename(Abundance_metabolism=n)
  ## Count the total number of genes per metabolism in each genome ##
  metabolism_counts_genome<-Kegg_long %>%
    select(all_of(met_features_y), Bin_name, KO, Abundance) %>%
    distinct() %>%
    filter(Abundance != "0") %>%
    group_by(Bin_name) %>%
    count(!!y_axis_features) %>%
    rename(Abundance_metabolism_genome=n)
  ############### Join and calculate the percentage #################
  Table_with_percentage<- left_join(
    metabolism_counts_genome, metabolism_counts_totals, 
    by=met_features_y) %>%
    mutate(Percentage = (Abundance_metabolism_genome * 100) / Abundance_metabolism ) %>%
    select(all_of(met_features_y), Bin_name, Percentage) %>%
    distinct() %>%
    drop_na() %>%
    left_join(other_data, by="Bin_name")
  ################################ Plot #############################
  Table_with_percentage_plot<-ggplot(Table_with_percentage,
                                     aes_string(x=x_axis,
                                                y=met_features_y,
                                                size=size_bubble,
                                                color=metadata_feature)) +
    geom_point(alpha=0.5) +
    scale_size(range = c(1,5)) +
    theme_linedraw() +
    theme(axis.text.x = element_text(size=6, 
                                     angle = 45, 
                                     hjust = 1, 
                                     vjust = 1),
        axis.text.y = element_text(size=5))

  return(Table_with_percentage_plot)
}


plot_bubble_percentage(kegg_map, 
                       metadata, 
                       "Bin_name",
                       Module,
                       "Percentage",
                       "Clades")

