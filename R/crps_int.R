#### Numerical integration ####
# Note: y can be a vector, all other inputs are scalars
crps_int <- function(y, pxxx, lower, upper, rel_tol = 1e-6){
  ind <- (y > upper) - (y < lower)
  out <- numeric(length(y))
  F1 <- function(x) pxxx(x)^2
  F2 <- function(x) (1-pxxx(x))^2
  if (any(ind == -1)) {
    out[ind == -1] <- sapply(which(ind == -1), function(i) {
      s1 <- lower - y[i]
      s2 <- integrate(F2, lower, upper, rel.tol = rel_tol)$value
      s1 + s2
    })
  } else if (any(ind == 0)) {
    out[ind == 0] <- sapply(which(ind == 0), function(i) {
      s1 <- integrate(F1, lower, y[i], rel.tol = rel_tol)$value
      s2 <- integrate(F2, y[i], upper, rel.tol = rel_tol)$value
      s1 + s2
    })
  } else if (any(ind == 1)) {
    out[ind == 1] <- sapply(which(ind == 1), function(i) {
      s1 <- integrate(F1, lower, upper, rel.tol = rel_tol)$value
      s2 <- y[i] - upper
      s1 + s2
    })
  }
  return(out)
}
