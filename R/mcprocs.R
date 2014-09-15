# Numerical integration for CRPS
crps.int <- function(F, y){
  s1 <- integrate(function(s) F(s)^2, -Inf, y)$value
  s2 <- integrate(function(s) (1-F(s))^2, y, Inf)$value
  s1 + s2
}

# CRPS of a (weighted) empirical distribution
crps.edf <- function(dat, y, w = NULL){
  
  ord <- order(dat)
  x <- dat[ord]
  n <- length(x)
  
  # Set uniform weights unless specified otherwise
  if (is.null(w)){
    dat <- .Internal(sort(dat,FALSE))
    out <- 2 / n^2 * sum((n * (y < dat) - 1:n + 0.5) * (dat - y))
  } 
  if (!is.null(w)) {
    if (length(w) != length(dat)) stop()
    w <- w[ord]
    p <- cumsum(w[1:(length(w)-1)])
  
  # Compute score
    if (y < x[1]){
    
      a <- x[2:n] - x[1:(n-1)]
      out <- sum(a*(p^2)) + x[1] - y
    
    } else if (y >= x[1] & y <= x[n]){
    
      ind <- min(which(x > y)) - 2
      a <- b <- rep(0, n-1)
      if (ind > 0){
        a[1:ind] <- x[2:(ind + 1)] - x[1:ind]
      }
      a[ind+1] <- y - x[ind+1]
      b[ind+1] <- x[ind+2] - y
      if (ind < (n-2)){
        b[(ind+2):(n-1)] <- x[(ind+3):(n)] - x[(ind+2):(n-1)]
      }
      out <- sum( a*(p^2) + b*((1-p)^2) )
    
    } else if (y > x[n]) {
    
      b <- x[2:n] - x[1:(n-1)]
      out <- sum(b*((1-p)^2)) + y - x[n]
    
    }
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
