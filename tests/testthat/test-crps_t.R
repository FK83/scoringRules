context("CRPS for t-distribution")

test_that("computed values are correct", {
  const <- 2.33875895
  
  expect_equal(crps_t(-3, 5), const)
  expect_equal(crps_t(-3 + .1, 5, location = .1), const)
  expect_equal(crps_t(-3 * .9, 5, scale = .9), const * .9)
  
  expect_equal(crps_tt(-3, 5), const)
  expect_equal(crps_tt(-3 + .1, 5, location = .1), const)
  expect_equal(crps_tt(-3 * .9, 5, scale = .9), const * .9)
  
  expect_equal(crps_ct(-3, 5), const)
  expect_equal(crps_ct(-3 + .1, 5, location = .1), const)
  expect_equal(crps_ct(-3 * .9, 5, scale = .9), const * .9)
  
  expect_equal(crps_gtct(-3, 5), const)
  expect_equal(crps_gtct(-3 + .1, 5, location = .1), const)
  expect_equal(crps_gtct(-3 * .9, 5, scale = .9), const * .9)
  
  const <- 0.566814455
  
  expect_equal(crps_tt(-1, 5, lower = -3, upper = 2), const)
  expect_equal(crps_gtct(-1, 5, lower = -3, upper = 2), const)
})
