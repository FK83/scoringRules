context("CRPS for normal distribution")

test_that("computed values are correct", {
  expect_equal(crps_norm(-3), 2.43657473)
  expect_equal(crps_tnorm(-3), 2.43657473)
  expect_equal(crps_cnorm(-3), 2.43657473)
  expect_equal(crps_gtcnorm(-3), 2.43657473)
  
  expect_equal(crps_norm(-3 + .1, location = .1), 2.43657473)
  expect_equal(crps_tnorm(-3 + .1, location = .1), 2.43657473)
  expect_equal(crps_cnorm(-3 + .1, location = .1), 2.43657473)
  expect_equal(crps_gtcnorm(-3 + .1, location = .1), 2.43657473)
  
  expect_equal(crps_norm(-3 * .9, scale = .9), 2.43657473 * .9)
  expect_equal(crps_tnorm(-3 * .9, scale = .9), 2.43657473 * .9)
  expect_equal(crps_cnorm(-3 * .9, scale = .9), 2.43657473 * .9)
  expect_equal(crps_gtcnorm(-3 * .9, scale = .9), 2.43657473 * .9)
})
