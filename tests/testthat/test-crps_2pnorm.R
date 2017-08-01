context("CRPS for two-piece-normal distribution")

test_that("computed values are correct", {
  const <- 2.80708311
  expect_equal(crps_2pnorm(-3, .7, 1.1), const)
  expect_equal(crps_2pnorm(-3 + .1, .7, 1.1, .1), const)
  
  const <- 2.17137978
  expect_equal(crps_2pnorm(3, .7, 1.1), const)
  expect_equal(crps_2pnorm(3 + .1, .7, 1.1, .1), const)
})
