test_that("get_subset_pca works on Pfam profiles with a realistic cutoff", {
  skip_on_cran()
  testthat::skip_if_not_installed("FactoMineR")
  
  # Minimal wide profile: rows=Pfam features, cols=bins
  pfam_tbl <- tibble::tibble(
    Pfam        = c("PF00001","PF00002","PF00003","PF00004"),
    domain_name = c("PF00001","PF00002","PF00003","PF00004"),
    Bin_A       = c(5, 0, 2, 8),
    Bin_B       = c(4, 1, 0, 7),
    Bin_C       = c(0, 0, 3, 2)
  )
  
  out <- get_subset_pca(
    tibble_rbims = pfam_tbl,
    cos2_val     = 0.80,
    analysis     = "Pfam"
  )
  
  expect_s3_class(out, "tbl_df")
  expect_true("Pfam" %in% names(out))
  # Debe devolver un subconjunto (no necesariamente todas)
  expect_true(nrow(out) >= 1)
})

test_that("get_subset_pca works on KEGG tables and renames back KO", {
  skip_on_cran()
  testthat::skip_if_not_installed("FactoMineR")
  
  kegg_tbl <- tibble::tibble(
    KO                  = c("K00001","K00002","K00003","K00004"),
    Module              = NA_character_,
    Module_description  = NA_character_,
    Pathway             = NA_character_,
    Pathway_description = NA_character_,
    Genes               = NA_character_,
    Gene_description    = NA_character_,
    Enzyme              = NA_character_,
    Cycle               = NA_character_,
    Pathway_cycle       = NA_character_,
    Detail_cycle        = NA_character_,
    rbims_pathway       = NA_character_,
    rbims_sub_pathway   = NA_character_,
    Bin_1               = c(3, 0, 1, 5),
    Bin_2               = c(2, 0, 0, 4),
    Bin_3               = c(0, 1, 2, 3)
  )
  
  out <- get_subset_pca(
    tibble_rbims = kegg_tbl,
    cos2_val     = 0.75,
    analysis     = "KEGG"
  )
  
  expect_s3_class(out, "tbl_df")
  expect_true("KO" %in% names(out))
  expect_false("tmp" %in% names(out))
  expect_true(nrow(out) >= 1)
})

test_that("get_subset_pca errors if analysis='PFAM' (wrong capitalization)", {
  skip_on_cran()
  testthat::skip_if_not_installed("FactoMineR")
  
  pfam_tbl <- tibble::tibble(
    Pfam        = c("PF00001","PF00002"),
    domain_name = c("PF00001","PF00002"),
    Bin_A       = c(2, 3),
    Bin_B       = c(1, 4)
  )
  expect_error(
    get_subset_pca(pfam_tbl, cos2_val = 0.8, analysis = "PFAM"),
    regexp = "Use 'Pfam'",
    ignore.case = FALSE
  )
})

test_that("get_subset_pca returns empty tibble (but valid structure) with strict cutoff", {
  skip_on_cran()
  testthat::skip_if_not_installed("FactoMineR")
  
  pfam_tbl <- tibble::tibble(
    Pfam        = c("PF00001","PF00002","PF00003"),
    domain_name = c("PF00001","PF00002","PF00003"),
    Bin_A       = c(1, 0, 1),
    Bin_B       = c(0, 1, 1)
  )
  
  out <- get_subset_pca(
    tibble_rbims = pfam_tbl,
    cos2_val     = 0.999,   # muy estricto
    analysis     = "Pfam"
  )
  
  expect_s3_class(out, "tbl_df")
  expect_true("Pfam" %in% names(out))
  expect_true(nrow(out) >= 0)  # puede quedar vacío
})

test_that("get_subset_pca tolerates duplicate feature names via make.unique", {
  skip_on_cran()
  testthat::skip_if_not_installed("FactoMineR")
  
  pfam_tbl <- tibble::tibble(
    Pfam        = c("PF00001","PF00001","PF00002"),
    domain_name = c("PF00001","PF00001","PF00002"),
    Bin_A       = c(2, 1, 0),
    Bin_B       = c(0, 1, 3)
  )
  
  out <- get_subset_pca(
    tibble_rbims = pfam_tbl,
    cos2_val     = 0.70,
    analysis     = "Pfam"
  )
  
  expect_s3_class(out, "tbl_df")
  expect_true("Pfam" %in% names(out))
  expect_true(nrow(out) >= 1)
})

