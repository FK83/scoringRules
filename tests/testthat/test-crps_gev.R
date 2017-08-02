context("CRPS for gev distribution")

FF <- function(shape) {
  if (shape == 0) {
    function(x) exp(-exp(-x))
  } else {
    function(x) {
      x <- 1 + shape * x
      x[x < 0] <- 0
      x <- x^(-1/shape)
      exp(-x)
    }
  }
}

test_that("computed values are correct", {
  const <- 0.276440963
  expect_equal(crps_gev(.3, 0), const)
  expect_equal(crps_gev(.3 + .1, 0, location = .1), const)
  expect_equal(crps_gev(.3 * .9, 0, scale = .9), const * .9)
  
  const <- 0.458044365
  expect_equal(crps_gev(.3, .7), const)
  expect_equal(crps_gev(.3 + .1, .7, location = .1), const)
  expect_equal(crps_gev(.3 * .9, .7, scale = .9), const * .9)
  
  const <- 0.207621488
  expect_equal(crps_gev(.3, -.7), const)
  expect_equal(crps_gev(.3 + .1, -.7, location = .1), const)
  expect_equal(crps_gev(.3 * .9, -.7, scale = .9), const * .9)
})
