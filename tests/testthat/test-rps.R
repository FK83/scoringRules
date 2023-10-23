context("Ranked probability score")

test_that("function works as expected", {
  expect_error(rps_probs(y = 5, p = c(.1, .5, .5)))
  expect_error(rps_probs(y = 3.1, p = c(.1, .5, .5)))
  expect_error(rps_probs(y = c(2, 3), p = c(.1, .5, .5)))
  K <- sample.int(1e3, size = 1)
  x <- runif(K)
  x <- x/sum(x)
  y <- sample.int(K, size = 1, prob = x)
  expect(rps_probs(y = y, x = x) > 0, 
         "RPS not strictly positive")
})
