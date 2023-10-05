test_that("EL STOP FUNCIONA", {
  expect_error(read_interpro(data_interpro = "inst/extdata/Interpro_test.tsv",
                             database="PFAM", profile = T), "PFAM va con minusculas!")
})
