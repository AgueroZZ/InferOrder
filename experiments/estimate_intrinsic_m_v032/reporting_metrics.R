# Use posterior occupancy consistently when reporting the effective dimension.
reported_effective_M <- function(result) {
  if (result$row$status != "success") return(NA_integer_)
  as.integer(sum(colMeans(result$selected$W) > design$effective_weight_tol))
}

reported_result_row <- function(result) {
  row <- result$row
  row$original_reported_M <- row$estimated_M
  row$selected_model_M <- if (row$status == "success") result$selected$M else NA_integer_
  row$effective_M <- reported_effective_M(result)
  row$estimated_M <- row$effective_M
  row
}
