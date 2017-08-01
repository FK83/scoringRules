context("CRPS for beta distribution")

test_that("computed values are correct", {
  const <- 0.0850102437
  expect_equal(crps_beta(.3, .7, 1.1), const)
  
  const <- 0.883206751
  expect_equal(crps_beta(-3, .7, 1.1, -5, 4), const)
})
