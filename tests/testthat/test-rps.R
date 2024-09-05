context("Ranked probability score")

test_that("function works as expected", {
  # Wrong input format
  expect_error(rps_probs(y = 5, x = c(.1, .5, .5)))
  expect_error(rps_probs(y = 3.1, x = c(.1, .5, .4)))
  expect_error(rps_probs(y = c(2, 3), x = c(.1, .5, .5)))
  # Correct inputs
  y <- 5
  x1 <- c(rep(0, y-1), 1)
  x2 <- rep(1/y, y)
  expect_equal(rps_probs(y = y, x = x1), 0)
  aux <- seq(from = 1/y, to = 1-1/y, by = 1/y)
  expect_equal(rps_probs(y = y, x = x2), sum(aux^2))
})
