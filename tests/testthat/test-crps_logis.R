context("CRPS for logistic distribution")

test_that("computed values are correct", {
  const <- 2.0971747
  
  expect_equal(crps_logis(-3), const)
  expect_equal(crps_logis(-3 + .1, location = .1), const)
  expect_equal(crps_logis(-3 * .9, scale = .9), const * .9)
  
  expect_equal(crps_tlogis(-3), const)
  expect_equal(crps_tlogis(-3 + .1, location = .1), const)
  expect_equal(crps_tlogis(-3 * .9, scale = .9), const * .9)
  
  expect_equal(crps_clogis(-3), const)
  expect_equal(crps_clogis(-3 + .1, location = .1), const)
  expect_equal(crps_clogis(-3 * .9, scale = .9), const * .9)
  
  expect_equal(crps_gtclogis(-3), const)
  expect_equal(crps_gtclogis(-3 + .1, location = .1), const)
  expect_equal(crps_gtclogis(-3 * .9, scale = .9), const * .9)
  
  const <- 0.513787837
  
  expect_equal(crps_tlogis(-1, lower = -3, upper = 2), const)
  expect_equal(crps_gtclogis(-1, lower = -3, upper = 2), const)
})
