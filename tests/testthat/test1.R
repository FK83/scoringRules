# Unit tests for crps.sample

expect_warning(crps.sample(runif(1000), 0.2, method = "kde", w = c(0.9, 0.1, rep(0, 998))))

# Check similarity of various forms to enter normal distribution

set.seed(22)
n <- 10
thresh <- 1e-4
y <- rnorm(n)
a <- matrix(0, n, 5)
b <- matrix(0, n, 4)

for (j in 1:n){

  a[j, 1] <- crps.parametric(family = "normal", parameters = list(m = 0, s = 1), y = y[j])
  a[j, 2] <- crps.parametric(family = "two-piece-normal", parameters = list(m = 0, s1 = 1, s2 = 1), y = y[j])
  a[j, 3] <- crps.parametric(family = "mixture-normal", parameters = list(m = 0, s = 1, w = 1), y = y[j])
  a[j, 4] <- crps.parametric(family = "t", parameters = list(m = 0, s = 1, df = 1e+7), y = y[j])
  expect_less_than(abs(a[j, 1]-a[j, 2]), thresh)
  expect_less_than(abs(a[j, 1]-a[j, 3]), thresh)
  expect_less_than(abs(a[j, 1]-a[j, 4]), thresh)
  
  b[j, 1] <- qs.parametric(family = "normal", parameters = list(m = 0, s = 1), y = y[j])
  b[j, 2] <- qs.parametric(family = "two-piece-normal", parameters = list(m = 0, s1 = 1, s2 = 1), y = y[j])
  b[j, 3] <- qs.parametric(family = "mixture-normal", parameters = list(m = 0, s = 1, w = 1), y = y[j])
  b[j, 4] <- qs.parametric(family = "t", parameters = list(m = 0, s = 1, df = 1e+5), y = y[j])
  expect_less_than(abs(b[j, 1]-b[j, 2]), thresh)
  expect_less_than(abs(b[j, 1]-b[j, 3]), thresh)
  expect_less_than(abs(b[j, 1]-b[j, 4]), thresh)
  
}
