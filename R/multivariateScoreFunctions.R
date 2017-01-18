################################################################################
###  multivariate scoring rules

################################################################################
# energy score

es_sample <- function(y, dat) {
  input <- list(y = y, dat = dat)
  check.multivsample(input)
  energyscoreC(y, dat)
}

################################################################################
# variogram score of order p

vs_sample <- function(y, dat, w = NULL,  p = 0.5) {
  input <- list(y = y, dat = dat)
  check.multivsample(input)
  d <- length(y)
  
  # additional input checks for weighting matrix w and order p
  if (is.null(w)) {
    w <- matrix(1, nrow = d, ncol = d)
  } else {
    if (!is.matrix(w)) {
      stop("'w' is not a matrix ")
    }
    if (any(dim(w) != d)) {
      stop("Dimensions of 'w' do not fit")
    }
    if (any(w < 0)) {
      stop("Weighting matrix 'w' contains negative values")
    }
  }

  if (!is.numeric(p) || length(p) != 1 ){
    stop("Order 'p' must be numeric of length 1")
  } else if (p < 0) {
    stop("Order 'p' must be positive")
  }
  
  
  out <- 0
  for (i in 1:d) {
    for (j in 1:d) {
      vdat <- mean(abs(dat[i, ] - dat[j, ])^p)
      vy <- abs(y[i] - y[j])^p
      out <- out + w[i, j] * (vy - vdat)^2 
    }
  }
  
  return(out)
}
