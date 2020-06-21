#' Reduce elements of genoprobs list by allelic series
#'
#' @param probs a list containing founder allele probabilities for a collection of markers
#' @param allelic_series a binary matrix of 0s and 1s indicating the allelic series
#' @return an object with class calc_genoprob.
#' @export
#' @details The output has two classes: list and calc_genoprob. Class 'calc_genoprob' is needed for input to qtl2::scan1, for example. Output length is same as that of inputted allele probs object.
reduce_probs <- function(probs, allelic_series){
  if (!is.list(probs)) stop("probs must be a list")
  if (!is.matrix(allelic_series) | !is.numeric(allelic_series)) stop("allelic_series must be a numeric matrix")
  if (ncol(allelic_series) < 2) stop("allelic_series must have at least 2 columns")
  out <- lapply(FUN = reduce_probs_array, X = probs, allelic_series = allelic_series)
  attr(out, "class") <- c("calc_genoprob", "list")
  return(out)
}

#' Reduce a genoprobs array by allelic series
#'
#' @param array a probs array (not a list)
#' @param allelic_series a matrix of 0s and 1s to indicate allelic series
#' @return a reduced array
#' @details The function multiplies every marker's probs matrix (from array inputted) by allelic_series and returns an array.
#' @export
reduce_probs_array <- function(array, allelic_series){
  out <- array(dim = c(dim(array)[[1]],
                       ncol(allelic_series),
                       dim(array)[[3]]))
  for (i in seq_along(1:(dim(array)[[3]]))){
    out[ , , i] <- array[ , , i] %*% allelic_series
  }
  rownames(out) <- rownames(array)
  colnames(out) <- paste0("allele", 1:ncol(out))
  dimnames(out)[[3]] <- dimnames(array)[[3]]
  return(out)
}
