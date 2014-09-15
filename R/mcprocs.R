# Numerical integration for CRPS
crps.int <- function(F, y){
  s1 <- integrate(function(s) F(s)^2, -Inf, y)$value
  s2 <- integrate(function(s) (1-F(s))^2, y, Inf)$value
  s1 + s2
}

# CRPS of a (weighted) empirical distribution
crps.edf <- function(dat, y, w = NULL){
  
  n <- length(dat)
  
  # Set uniform weights unless specified otherwise
  if (is.null(w)){
    x <- .Internal(sort(dat, FALSE))
    out <- 2 / n^2 * sum((n * (y < x) - 1:n + 0.5) * (x - y))
  } else {
    if (length(w) != n) stop()
    ord <- order(dat)
    x <- dat[ord]
    w <- w[ord]
    p <- c(0, cumsum(w))
    out <- 2 * sum((w * (y < x) - 0.5 * (p[2:(n+1)]^2 - p[1:n]^2)) * (x - y))
  }
  
  return(out)
}

# Estimate the CRPS of a sample of data (based on kernel representation)
crps.krep <- function(dat, y){
  
  n <- length(dat)
  shift <- as.integer((n-1)/2)
  shift <- c((shift+1):n, 1:shift)
  crps <- mean(abs(y-dat)) - 0.5*mean(abs(dat- dat[shift]))
  return(crps)
  
}

# CRPS for a mixture of normals
crps.mixn <- function(m, v, y, w = NULL){
  n <- length(m)
  if (is.null(w)) w <- rep(1/n, n)
  -crpsmixnC(w, m, sqrt(v), y)
}

# CRPS for a t distribution
crps.t <- function(m, v, df, y){
  s <- sqrt(v)
  z <- (y-m)/s
  c1 <- (2*s*dt(z, df)) * ((df + z^2)/(df - 1)) + z*s*(2*pt(z, df) - 1)
  c2 <- ( s*(4*sqrt(df))/(df - 1) ) * (beta(0.5, df - 0.5)/(beta(0.5, 0.5*df)^2))
  c1 - 0.5*c2
}

# CRPS via kernel density estimation
crps.kdens <- function(dat, y, bw = NULL){
  n <- length(dat)
  if (is.null(bw)){
    v <- rep(bw.SJ(dat)^2, n)
  } else {
    v <- rep(bw^2, n)
  }
  return(crps.mixn(dat, v, y))
}
