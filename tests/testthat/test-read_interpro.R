test_that("EL STOP FUNCIONA", {
  path <- system.file("extdata", "Interpro_test.tsv", package = "rbims")
  expect_error(
    read_interpro(data_interpro = path, database = "PFAM", profile = TRUE),
    regexp = "P letter must be in capital followed by lower-case letters e.g 'Pfam'"
  )
})
