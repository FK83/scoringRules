context("CRPS for mixture of normal distributions")

test_that("computed values are correct", {
  expect_equal(crps_mixnorm(-1,
                            matrix(-1.4, 1, 1),
                            matrix(.9, 1, 1),
                            matrix(1, 1, 1)),
               crps_norm(-1,-1.4, .9))
  expect_equal(crps_mixnorm_int(-1,
                                matrix(-1.4, 1, 1),
                                matrix(.9, 1, 1),
                                matrix(1, 1, 1),
                                rel_tol = 1e-10),
               crps_norm(-1,-1.4, .9))
})
