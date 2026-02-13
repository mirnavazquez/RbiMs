#' Calculate pathway-level directional bias (binomial test)
#'
#' Aggregates KO-level effects into pathway-level directionality statistics.
#' This is useful when many genes in a pathway consistently shift in the same
#' direction even if individual KO tests are not significant.
#'
#' Effect direction convention:
#' - `effect > 0` means enrichment toward the "positive" group in your contrast
#'   (depends on how the contrast is built in the discriminant pipeline).
#'
#' @param disc_obj The discriminant object (attr(, "rbims_disc")) containing
#'   at least `consensus` with columns `feature` and `effect`.
#' @param ko_dictionary A KO dictionary as produced by `make_ko_dictionary()`,
#'   with columns `KO` and `rbims_pathway`.
#' @param pathways Character vector of pathways to test (values from `rbims_pathway`).
#' @param p_null Null proportion for directionality. Default 0.5.
#' @param p_alternative Alternative hypothesis for directionality p-value.
#'   Default "greater" tests whether prop(effect>0) > 0.5.
#' @param ci_two_sided Logical. If TRUE, compute a two-sided 95% exact CI
#'   for the proportion (Clopper–Pearson), while keeping a directional p-value
#'   defined by `p_alternative`. This is often preferable for figures.
#' @param conf_level Confidence level for CI. Default 0.95.
#' @param drop_na_effect If TRUE, drop rows with NA effect before calculations.
#'
#' @return A tibble with `rbims_pathway`, `x_pos`, `n`, `prop_pos`,
#'   `ci_low`, `ci_high`, `p_value`, `p_adj_fdr`.
#' @export
calc_pathway_directional_bias <- function(disc_obj,
                                          ko_dictionary,
                                          pathways = NULL,
                                          p_null = 0.5,
                                          p_alternative = c("greater", "less", "two.sided"),
                                          ci_two_sided = TRUE,
                                          conf_level = 0.95,
                                          drop_na_effect = TRUE) {
  
  p_alternative <- match.arg(p_alternative)
  
  stopifnot(is.list(disc_obj))
  stopifnot(!is.null(disc_obj$consensus))
  stopifnot(is.data.frame(disc_obj$consensus))
  stopifnot(all(c("feature", "effect") %in% names(disc_obj$consensus)))
  
  stopifnot(is.data.frame(ko_dictionary))
  stopifnot(all(c("KO", "rbims_pathway") %in% names(ko_dictionary)))
  
  # Annotate consensus with pathway labels
  cons_annot <- disc_obj$consensus %>%
    dplyr::left_join(ko_dictionary, by = c("feature" = "KO"))
  
  if (drop_na_effect) {
    cons_annot <- cons_annot %>% dplyr::filter(!is.na(effect))
  }
  
  # Keep only rows with a pathway
  cons_annot <- cons_annot %>% dplyr::filter(!is.na(rbims_pathway))
  
  # If user did not provide pathways, use all observed pathways
  if (is.null(pathways)) {
    pathways <- sort(unique(cons_annot$rbims_pathway))
  }
  
  cons_annot <- cons_annot %>% dplyr::filter(rbims_pathway %in% pathways)
  
  res <- lapply(pathways, function(pw) {
    
    tmp <- cons_annot %>% dplyr::filter(rbims_pathway == pw)
    
    x <- sum(tmp$effect > 0, na.rm = TRUE)
    n <- nrow(tmp)
    
    if (n == 0) {
      return(tibble::tibble(
        rbims_pathway = pw,
        x_pos = NA_integer_,
        n = 0L,
        prop_pos = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_,
        p_value = NA_real_
      ))
    }
    
    # CI: optionally two-sided (more standard for plots)
    if (isTRUE(ci_two_sided)) {
      bt_ci <- stats::binom.test(x = x, n = n, conf.level = conf_level, alternative = "two.sided")
      ci_low <- bt_ci$conf.int[1]
      ci_high <- bt_ci$conf.int[2]
    } else {
      # CI aligned with p_alternative (often ends at 1.0 for "greater")
      bt_ci <- stats::binom.test(x = x, n = n, conf.level = conf_level, alternative = p_alternative)
      ci_low <- bt_ci$conf.int[1]
      ci_high <- bt_ci$conf.int[2]
    }
    
    # Directional p-value
    bt_p <- stats::binom.test(x = x, n = n, p = p_null, alternative = p_alternative)
    
    tibble::tibble(
      rbims_pathway = pw,
      x_pos = x,
      n = n,
      prop_pos = x / n,
      ci_low = ci_low,
      ci_high = ci_high,
      p_value = bt_p$p.value
    )
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(p_adj_fdr = stats::p.adjust(.data$p_value, method = "fdr")) %>%
    dplyr::arrange(.data$p_value)
  
  res
}