#' Simulate one trait with independent errors
#'
#' @param aprobs founder allele probabilities matrix
#' @param allelic_series a binary matrix of 0s and 1s indicating the allelic series
#' @param allelic_effects a vector of effects for the alleles
#' @param error_variance error variance for the traits
#' @export
#' @return a vector of values for a single trait

sim1 <- function(aprobs, allelic_series, allelic_effects, error_variance = 1){
  return(aprobs %*% allelic_series %*% allelic_effects + stats::rnorm(n = nrow(aprobs), sd = sqrt(error_variance)))
}
