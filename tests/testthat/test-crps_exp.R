context("CRPS for exponential distribution")

test_that("computed values are correct", {
  const <- 1.20701837
  expect_equal(crps_exp(3, .7), const)
  expect_equal(crps_expM(3, scale = 1 / 0.7), const)
  
  const <- 2.11137878
  expect_equal(crps_expM(3, scale = 1 / 0.7, mass = .6), const)
})
