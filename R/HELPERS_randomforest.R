#' @keywords internal
rf_importance_features <- function(mat_prop, group){
  rf <- randomForest::randomForest(
    x = t(mat_prop),
    y = as.factor(group),
    importance = TRUE
  )
  tibble::tibble(
    feature = rownames(rf$importance),
    rf_importance = rf$importance[, 1]
  )
}
