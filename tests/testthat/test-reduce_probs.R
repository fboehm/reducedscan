library(testthat)
library(reducedscan)

# use real qtl2 outputs to test reduce_probs

# read data
is_inst <- function(pkg) {
  nzchar(system.file(package = pkg))
}
qtl2_indic <- is_inst("qtl2")

ghurl <- "https://raw.githubusercontent.com/rqtl/qtl2data/master/DOex/DOex_genoprobs.rds"
probs <- readRDS(gzcon(url(ghurl)))
aprobs <- qtl2::genoprob_to_alleleprob(probs)
# define allelic_series
  allelic_series <- matrix(c(1, 0, 0,
                             1, 0, 0,
                             0, 1, 0,
                             0, 1, 0,
                             0, 1, 0,
                             0, 0, 1,
                             0, 0, 1,
                             0, 0, 1 ),
                           nrow = 8, byrow = TRUE
  )

# tests

test_that("reduce_probs returns a list of proper length", {
  expect_length(reduce_probs(probs = aprobs, allelic_series = allelic_series), 3)
  expect_equal(class(reduce_probs(probs = aprobs, allelic_series = allelic_series)), c("calc_genoprob", "list"))
})

test_that("reduce_probs returns list containing arrays of correct dimension", {
  expect_equal(dim(reduce_probs(probs = aprobs, allelic_series = allelic_series)[[1]]), c(dim(aprobs[[1]])[[1]], ncol(allelic_series), dim(aprobs[[1]])[[3]]))
})

