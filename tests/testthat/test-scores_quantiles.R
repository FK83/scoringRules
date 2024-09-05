context("Quantile and interval scores")

test_that("computed values are correct", {
  expect_equal(qs_quantiles(y = 1, x = 0, alpha = .1), .1)
  expect_equal(qs_quantiles(y = 1, x = 0, alpha = .4), 
               4*qs_quantiles(y = 1, x = 0, alpha = .1))
  expect_equal(qs_quantiles(y = -1, x = 0, alpha = .1), .9)

  expect_equal(ints_quantiles(y = 1, x_lower = -1, x_upper = 1, target_coverage = .8), 
               2)
  expect_equal(ints_quantiles(y = 1, x_lower = -1, x_upper = 1, target_coverage = .8), 
               ints_quantiles(y = 1, x_lower = -1, x_upper = 1, target_coverage = .5))
  expect_equal(ints_quantiles(y = 2, x_lower = -1, x_upper = 1, target_coverage = .8), 
               12)
  expect_equal(ints_quantiles(y = 2, x_lower = -1, x_upper = 1, target_coverage = .5), 
               6)
  
  dat <- qnorm(seq(from = .01, to = .99, by = .01))
  expect_equal(ints_sample(y = 1, dat = dat, target_coverage = .6), 
               5*(qs_sample(y = 1, dat = dat, alpha = .2) + 
                    qs_sample(y = 1, dat = dat, alpha = .8)))
  
  expect_error(ints_sample(y = 1, dat = dat, target_coverage = 1.1))
  expect_error(ints_sample(y = 1, dat = dat, target_coverage = .8, 
                           type = "a"))
  
})
