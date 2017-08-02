context("CRPS for log-normal distribution")

FF <- function(x) plnorm(x, .1, .9)

test_that("computed values are correct", {
  const <- 1.13552645
  expect_equal(crps_lnorm(3, .1, .9), const)
})
