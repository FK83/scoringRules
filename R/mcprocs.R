# Simulate Data from the Fox/West Model
sim.fw <- function(f, s, n, ndraws){
  # Compute parameters of marginal distribution (inverse Gamma)
  beta.ig <- c(0.5*(n+2), 0.5*n*s)
  # Compute parameters of corresponding t distribution
  beta.t <- c(beta.ig[2]/beta.ig[1], 2*beta.ig[1]) # variance, degrees of freedom
  # Simulate data
  sig2 <- rep(0, ndraws)
  sig2[1] <- 1
  psi <- 1/rgamma(ndraws, 0.5*(n+3), 0.5*n*s*(1-f^2))
  v <- f + sqrt( psi/(n*s) ) * rnorm(ndraws)
  for (j in 2:ndraws){
    sig2[j] <- psi[j] + (v[j]^2) * sig2[j-1]
  }  
  # Return list
  return(list(par.ig = beta.ig, sim = sig2))
}

# Make the ergodic dist for our simulation design
ergdist <- function(s, n){
  # Make t distribution parameters
  v <- n*s/(n + 2)
  df <- n + 2
  # Compute properties of ergodic dist
  pdf <- function(z){
    sapply(z, function(l) dt(l/sqrt(v), df)/sqrt(v))
  }  
  cdf <- function(z){
    sapply(z, function(l) pt(l/sqrt(v), df))
  } 
  list(pdf = pdf, cdf = cdf, v = v, df = df)
}

# Numerical integration for CRPS
crps.int <- function(F, y, subdivisions = 1000L, stop.on.error = TRUE){
  s1 <- integrate(function(s) F(s)^2, -Inf, y, subdivisions = subdivisions, stop.on.error = stop.on.error)$value
  s2 <- integrate(function(s) (1-F(s))^2, y, Inf, subdivisions = subdivisions, stop.on.error = stop.on.error)$value
  s1 + s2
}

# CRPS of a (weighted) empirical distribution
crps.edf <- function(dat, y, w = NULL){
  n <- length(dat)
  # Set uniform weights unless specified otherwise
  if (is.null(w)){
    x <- .Internal(sort(dat, FALSE))
    out <- sapply(y, function(s) 2 / n^2 * sum((n * (s < x) - 1:n + 0.5) * (x - s)))
  } else {
    if (length(w) != n) stop()
    ord <- order(dat)
    x <- dat[ord]
    w <- w[ord]
    p <- c(0, cumsum(w))
    out <- sapply(y, function(s) 2 * sum((w * (s < x) - 0.5 * (p[2:(n+1)]^2 - p[1:n]^2)) * (x - s)))
  }
  return(out)
}

# Estimate the CRPS of a sample of data (based on kernel representation)
crps.krep <- function(dat, y){
  
  n <- length(dat)
  shift <- (n-1) %/% 2
  shift <- c((shift+1):n, 1:shift)
  crps <- mean(abs(y-dat)) - 0.5*mean(abs(dat- dat[shift]))
  return(crps)
  
}

# CRPS for a mixture of normals
crps.mixn = function (m, v, y, w = NULL, exact = TRUE, subdivisions = 1000L, stop.on.error = TRUE){
  n <- length(m)
  if (is.null(w)) 
    w <- rep(1/n, n)
  if (exact == TRUE){
    return(-crpsmixnC(w, m, sqrt(v), y))
  } else {
    Fmix = function(z){
      sapply(z, function(r) sum(w*pnorm((r-m)/sqrt(v))))
    }
    return(crps.int(Fmix, y, subdivisions = subdivisions, stop.on.error = stop.on.error))
  }
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
crps.kdens = function(dat, y, bw = NULL, exact = TRUE, subdivisions = 1000L, stop.on.error = TRUE){
  n <- length(dat)
  if (is.null(bw)) {
    v <- rep(bw.SJ(dat)^2, n)
  }
  else {
    v <- rep(bw^2, n)
  }
  return(crps.mixn(dat, v, y, exact = exact, subdivisions = subdivisions, stop.on.error = stop.on.error))
}

# Log Score for Mixture of Normals
ls.mixn <- function(m, v, y, w = NULL) {
  n <- length(m)
  if (is.null(w)){
    w <- rep(1/n, n)
  }
  lsmixnC(m, sqrt(v), w, y)
}

# Log Score for Kernel Density Estimate
ls.kdens <- function(dat, y, bw = NULL) {
  if (is.null(bw)) 
    bw <- bw.SJ(dat)
  ls.mixn(dat, rep(bw^2, length(dat)), y)
}