context("CRPS for two-piece-exponential distribution")

test_that("computed values are correct", {
  const <- 2.72138251
  expect_equal(crps_2pexp(-3, .7, 1.1), const)
  
  const <- 2.00181206
  expect_equal(crps_2pexp(3, .7, 1.1), const)
})
