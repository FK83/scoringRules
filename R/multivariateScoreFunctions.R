################################################################################
###  multivariate scoring rules

################################################################################
# energy score

euclnorm <- function(x) {
  sqrt(sum(x %*% x))
}

es_sample <- function(y, dat, use_dist = FALSE) {
  input <- list(y = y, dat = dat)
  check.multivsample(input)
  if (length(use_dist) > 1 | !is.logical(use_dist)) {
    stop("'use_dist' must be logical of length 1")
  }
  
  out <- numeric()
  m <- dim(dat)[2]
  s1 <- sum(sapply(1:m, function(i){euclnorm(dat[, i] - y) } ))
  if (use_dist) {
    s2 <- 2*sum(dist(t(dat)))
  } else {
    s2 <- 0
    for (j in 1:m) {
      for (k in j:m) {
        s2 <- s2 +2 * euclnorm(dat[, j] - dat[, k]) 
      }
    }
  }
  out <- (s1 / m) - s2 / (2 * m * m)
  
  return(out)
}

################################################################################
### variogram score of order p

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