library(reducedscan)
library(testthat)

# pre-test code

## make array

probs <- array(data = runif(n = 800, max = 1 / 8), dim = c(20, 8, 5))


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
test_that("output is array with correct dimensions", {
  expect_equal(dim(reduce_probs_array(probs, allelic_series)),
               c(dim(probs)[[1]], ncol(allelic_series), dim(probs)[[3]])
  )
})

