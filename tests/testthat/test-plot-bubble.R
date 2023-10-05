test_that("EL STOP FUNCIONA", {
  expect_error(plot_bubble(interpro_pfam_profile, 
                                       y_axis=Pfam, 
                                       x_axis=Bin_name, 
                                       analysis = "INTERPRO", 
                                       data_experiment = metadata, 
                                       color_character = Clades), 
               "calc must have a value between Abundance, Binary or Percentage")
})
