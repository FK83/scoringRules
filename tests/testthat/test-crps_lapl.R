context("CRPS for Laplace distribution")

FF <- function(x) {
  0.5 + 0.5 * sign(x) * pexp(abs(x))
}

test_that("computed values are correct", {
  const <- 2.29978707
  expect_equal(crps_lapl(-3), const)
  expect_equal(crps_lapl(-3 + .1, location = .1), const)
  expect_equal(crps_lapl(-3 * .9, scale = .9), .9 * const)
})
