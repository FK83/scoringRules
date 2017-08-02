context("CRPS for log-Laplace distribution")

FF <- function(x) {
  0.5 + 0.5 * sign(log(x) - .1) * pexp(abs(log(x) - .1) / .9)
}

test_that("computed values are correct", {
  const <- 1.16202051
  expect_equal(crps_llapl(3, .1, .9), const)
})
