test_that("calc_abundance cuenta ocurrencias por Bin y feature (tmp) y respeta col_rename", {
  skip_on_cran()
  ds <- tibble::tibble(
    Bin_name = c("Bin_01_scaffold1","Bin_01_scaffold2","Bin_02_scaffold3"),
    KO       = c("K1","K1","K2")
  )
  # Nota: analysis="KEGG" usa columna 'KO' como feature (tmp)
  out1 <- calc_abundance(ds, analysis = "KEGG")
  expect_true(all(c("Bin_name","Scaffold_name","tmp","Abundance") %in% names(out1)))
  # Abundance por Bin y tmp:
  # Bin_01 tiene K1 dos veces, Bin_02 tiene K2 una vez
  a_01 <- out1$Abundance[out1$Bin_name=="Bin_01" & out1$tmp=="K1"][1]
  a_02 <- out1$Abundance[out1$Bin_name=="Bin_02" & out1$tmp=="K2"][1]
  expect_equal(a_01, 2)
  expect_equal(a_02, 1)
  
  # col_rename cambia nombre de 'tmp'
  out2 <- calc_abundance(ds, analysis = "KEGG", col_rename = "KO_id")
  expect_true("KO_id" %in% names(out2))
  expect_false("tmp" %in% names(out2))
})
