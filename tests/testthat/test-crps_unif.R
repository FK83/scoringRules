context("CRPS for uniform distribution")

FF <- function(min, max, lmass, umass) {
  function(x) lmass + (1 - (lmass + umass)) * punif(x, min ,max)
}

test_that("computed values are correct", {
  const <- 0.123333333
  expect_equal(crps_unif(.3), const)
  
  const <- 0.248333333
  expect_equal(crps_unif(.3, 0, 1, .1, .4), const)
  
  const <- 1.01666667
  expect_equal(crps_unif(-1, -3, 2, .1, .4), const)
})
