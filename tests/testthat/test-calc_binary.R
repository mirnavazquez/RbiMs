test_that("calc_binary en modo binario devuelve presencia/ausencia por vía y MAG", {
  skip_on_cran()
  tib <- tibble::tibble(
    rbims_pathway = c("A","A","B","B"),
    KO            = c("K1","K2","K3","K4"),
    Bin_10        = c(1,0,0,0),
    Bin_12        = c(0,1,1,0)
  )
  out <- calc_binary(tib, rbims_pathway, binary = TRUE)
  expect_true(all(c("rbims_pathway","Bin_name","Presence_absence") %in% names(out)))
  
  getB <- function(bin, v) out$Presence_absence[out$Bin_name==bin & out$rbims_pathway==v]
  expect_equal(getB("Bin_10","A"), 1L)
  expect_equal(getB("Bin_10","B"), 0L)
  expect_equal(getB("Bin_12","A"), 1L)
})

test_that("calc_binary en modo abundancia suma Abundance por vía y MAG", {
  skip_on_cran()
  tib <- tibble::tibble(
    rbims_pathway = c("A","A","A"),
    KO            = c("K1","K2","K3"),
    Bin_10        = c(2,0,3),   # suma 5
    Bin_12        = c(0,1,0)    # suma 1
  )
  out <- calc_binary(tib, rbims_pathway, binary = FALSE)
  expect_true(all(c("rbims_pathway","Bin_name","Abundance") %in% names(out)))
  getA <- function(bin) out$Abundance[out$Bin_name==bin & out$rbims_pathway=="A"]
  expect_equal(getA("Bin_10"), 5)
  expect_equal(getA("Bin_12"), 1)
})

test_that("calc_binary puede unir metadata y devolver tabla extendida con metabolism=TRUE", {
  skip_on_cran()
  tib <- tibble::tibble(
    rbims_pathway = c("A","B"),
    KO            = c("K1","K2"),
    Bin_01        = c(1,0),
    Bin_02        = c(0,1)
  )
  meta <- tibble::tibble(Bin_name = c("Bin_01","Bin_02"), Class=c("Alpha","Gamma"))
  out <- calc_binary(tib, rbims_pathway, data_experiment = meta, binary = TRUE, metabolism = TRUE)
  expect_true("Class" %in% names(out)) # metadatos unidos
  expect_true(all(c("KO","Bin_name") %in% names(out))) # metabolism=TRUE reanexa info
})
