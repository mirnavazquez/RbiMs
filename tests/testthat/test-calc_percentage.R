test_that("calc_percentage calcula x/k*100 por MAG y vía", {
  skip_on_cran()
  # Vía A tiene KOs: K1,K2; vía B: K3,K4
  # Bin_10 tiene ambos KOs de A -> 100% en A, ninguno de B -> 0 (NA en plot)
  # Bin_12 tiene 1/2 de A -> 50%, 2/2 de B -> 100%
  tib <- tibble::tibble(
    rbims_pathway = c("A","A","B","B"),
    KO            = c("K1","K2","K3","K4"),
    Bin_10        = c(1,1,0,0),
    Bin_12        = c(1,0,1,1)
  )
  
  out <- calc_percentage(tib, rbims_pathway)
  
  # columnas esperadas
  expect_true(all(c("Pathway_number_of_total_elements","rbims_pathway","Bin_name","Percentage") %in% names(out)))
  
  # valores por MAG & vía (redondeamos solo para el assert de igualdad)
  getP <- function(bin, v) out$Percentage[out$Bin_name==bin & out$rbims_pathway==v]
  expect_equal(round(getP("Bin_10","A")), 100)
  expect_equal(round(getP("Bin_12","A")), 50)
  expect_equal(round(getP("Bin_12","B")), 100)
  
  # ceros no cuentan
  expect_true(!any(out$Percentage[out$Bin_name=="Bin_10" & out$rbims_pathway=="B"] > 0, na.rm = TRUE))
})

test_that("calc_percentage ignora Abundance==0 y mantiene decimales", {
  skip_on_cran()
  tib <- tibble::tibble(
    rbims_pathway = c("C","C","C","C","C"),
    KO            = c("K1","K2","K3","K4","K5"),
    Bin_01        = c(1,1,0,0,0),   # 2/5 => 40.0
    Bin_02        = c(1,0,1,0,0)    # 2/5 => 40.0 (mismo)
  )
  out <- calc_percentage(tib, rbims_pathway)
  vals <- out$Percentage[out$Bin_name %in% c("Bin_01","Bin_02")]
  expect_true(all(abs(vals - 40.0) < 1e-8))
})
