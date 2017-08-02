context("CRPS for gamma distribution")

test_that("computed values are correct", {
  const <- 0.399009355
  expect_equal(crps_gamma(.2, 1.1), const)
  expect_equal(crps_gamma(.2 * .9, 1.1, scale = .9), const * .9)
})
