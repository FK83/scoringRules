context("CRPS for log-logistic distribution")

FF <- function(x) plogis(log(x), .1, .9)

test_that("computed values are correct", {
  const <- 1.13295277
  expect_equal(crps_llogis(3, .1, .9), const)
})
