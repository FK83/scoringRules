context("CRPS for gpd distribution")

FF <- function(shape) {
  if (shape == 0) {
    function(x) {
      x <- exp(-x)
      x[x > 1] <- 1
      1 - x
    }
  } else {
    function(x) {
      x <- 1 + shape * x
      x[x < 0] <- 0
      x <- x^(-1/shape)
      x[x > 1] <- 1
      1 - x
    }
  }
}

test_that("computed values are correct", {
  const <- 0.281636441
  expect_equal(crps_gpd(.3, 0), const)
  expect_equal(crps_gpd(.3 + .1, 0, location = .1), const)
  expect_equal(crps_gpd(.3 * .9, 0, scale = .9), const * .9)
  
  const <- 0.546254141
  expect_equal(crps_gpd(.3, .7), const)
  expect_equal(crps_gpd(.3 + .1, .7, location = .1), const)
  expect_equal(crps_gpd(.3 * .9, .7, scale = .9), const * .9)
  
  const <- 0.157583485
  expect_equal(crps_gpd(.3, -.7), const)
  expect_equal(crps_gpd(.3 + .1, -.7, location = .1), const)
  expect_equal(crps_gpd(.3 * .9, -.7, scale = .9), const * .9)
})
