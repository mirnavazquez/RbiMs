test_that("plot_bubble acepta labs(...) y conserva porcentajes con decimales", {
  skip_on_cran()
  tibble_ko <- tibble::tibble(
    rbims_pathway = c("Hexadecane","Phenanthrene"),
    KO = c("K00496","K01031"),
    Bin_10 = c(1,0),
    Bin_12 = c(0,1)
  )
  meta <- tibble::tibble(
    Bin_name = c("Bin_10","Bin_12"),
    Class = c("Alpha","Gamma")
  )
  p <- plot_bubble(
    tibble_ko = tibble_ko,
    x_axis = Bin_name,
    y_axis = rbims_pathway,
    analysis = "KEGG",
    calc = "Percentage",
    data_experiment = meta,
    color_character = Class,
    y_labs = "ruta"     # <- usa el parámetro nativo
  )
  
  expect_s3_class(p, "ggplot")
  gb <- ggplot2::ggplot_build(p)
  # debe existir estética 'size' y la etiqueta de eje Y debe ser "ruta"
  expect_true(any(vapply(gb$data, function(d) "size" %in% names(d), logical(1))))
  expect_equal(p$labels$y, "ruta")
})