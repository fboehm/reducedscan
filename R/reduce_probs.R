#' Reduce elements of genoprobs list by allelic series
#'
#'
#'  @param genoprobs a list containing founder allele probabilities for a collection of markers
#'  @param allelic_series a binary matrix of 0s and 1s indicating the allelic series
reduce_probs <- function(probs, allelic_series){
  if (!is.list(probs)) stop("probs must be a list")
  if (!is.matrix(allelic_series) | !is.numeric(allelic_series)) stop("allelic_series must be a numeric matrix")
  if (ncol(allelic_series) < 2) stop("allelic_series must have at least 2 columns")
  out <- lapply(FUN = reduce_probs_array, X = probs)
  attr(out, "class", exact = TRUE) <- c(attr(out), attr(probs))
  return(out)
}

#' Reduce a genoprobs array by allelic series
#'
#' @param array a probs array (not a list)
#' @param allelic_series a matrix of 0s and 1s to indicate allelic series
#' @return a reduced array
#' @details The function multiplies every marker's probs matrix (from array inputted) by allelic_series and returns an array.
reduce_probs_array <- function(array, allelic_series){
  out <- array(dim = c(dim(array)[[1]], ncol(allelic_series), dim(array)[[3]]))
  for (i in seq_along(1:(dim(array)[[3]]))){
    out[ , , i] <- array[ , , i] %*% allelic_series
  }
  return(out)
}
