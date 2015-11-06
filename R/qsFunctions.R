################################################################################
### bounded interval

# uniform
qs.unif <- function(y, min, max) 2*dunif(y, min, max) - 1/(max-min)

# beta
qs.beta <- function(y, shape1, shape2) {
  c1 <- 2 * dbeta(y, shape1, shape2)
  c2 <- beta(2*shape1 - 1, 2*shape2 - 1) / beta(shape1, shape2)^2
  return(c1 - c2)
}


################################################################################
### real line

# laplace
qs.lapl <- function(y, location, scale) {
  c1 <- 2*flapl(y, location, scale)
  c2 <- 1/(4*scale)
  return(c1 - c2)
}

# logistic
qs.logis <- function(y, location, scale) {
  c1 <- 2*dlogis(y, location, scale)
  c2 <- 1/(6*scale)
  return(c1 - c2)
}

# normal
qs.norm <- function(y, mean, sd) {
  c1 <- 2*dnorm(y, mean, sd)
  c2 <- 1/(2*sd*sqrt(pi))
  return(c1 - c2)
}

# mixture of normals
qs.mixnorm <- function(y, m, s, w) sapply(seq_along(y), function(i) qsmixnC(w[i, ], m[i, ], s[i, ], y[i]))

# two-piece-normal
qs.2pnorm <- function(y, m, s1, s2) {
  c1 <- 2*f2pnorm(y, m, s1, s2)
  c2 <- (s1 + s2) / (sqrt(pi) * (s1 + s2)^2)
  return(c1 - c2)
}

# t
qs.t <- function(y, df, location, scale) {
  c1 <- 2*ft(y, df, location, scale)
#  c2 <- 1/df^(1.5) / sqrt(pi) / scale * gamma(df + 0.5) / gamma(df) * gamma((df + 1)/2)^2 / gamma(df/2)^2
#  c3 <- 2^(3-4*df) * pi / scale / df^(1.5) * gamma(2*df) / gamma(df/2)^4
  c2 <- 2^(3-4*df) * pi / scale / df^(1.5) / beta(df, df) / beta(df/2, df/2)^2
  ind <- !is.finite(c2)
  if (any(ind)) {
    c2[ind] <- rep(1/(2*scale*sqrt(pi)), len = length(c2))[ind]
  }
  
  return(c1 - c2)
}

################################################################################
### non-negative

# exponential
qs.exp <- function(y, rate) 2*dexp(y, rate) - 0.5*rate

# gamma
qs.gamma <- function(y, shape, scale) {
  c1 <- 2*dgamma(y, shape, scale = scale)
  c2 <- 2^(1-2*shape) / scale / (2*shape - 1) / beta(shape, shape)
  return(c1 - c2)
}

# log-laplace
qs.llapl <- function(y, locationlog, scalelog) {
  warning("No closed form expression implemented - numerical integration used.")
  n <- max(lengths(list(y, locationlog, scalelog)))
  y <- rep(y, len = n)
  locationlog <- rep(locationlog, len = n)
  scalelog <- rep(scalelog, len = n)
  
  c1 <- 2*fllapl(y, locationlog, scalelog)
  c2 <- sapply(1:n, function(i) {
    d <- function(x) fllapl(x, locationlog[i], scalelog[i])^2
    integrate(d, 0, Inf)$value
  })
  
  return(c1 - c2)
}

# log-logistic
qs.llogis <- function(y, locationlog, scalelog) {
  warning("No closed form expression implemented - numerical integration used.")
  n <- max(lengths(list(y, locationlog, scalelog)))
  y <- rep(y, len = n)
  locationlog <- rep(locationlog, len = n)
  scalelog <- rep(scalelog, len = n)
  
  c1 <- 2*fllogis(y, locationlog, scalelog)
  c2 <- sapply(1:n, function(i) {
    d <- function(x) fllogis(x, locationlog[i], scalelog[i])^2
    integrate(d, 0, Inf)$value
  })
  
  return(c1 - c2)
}

# log-normal
qs.lnorm <- function(y, meanlog, sdlog) {
  warning("No closed form expression implemented - numerical integration used.")
  n <- max(lengths(list(y, meanlog, sdlog)))
  y <- rep(y, len = n)
  meanlog <- rep(meanlog, len = n)
  sdlog <- rep(sdlog, len = n)
  
  c1 <- 2*dlnorm(y, meanlog, sdlog)
  c2 <- sapply(1:n, function(i) {
    d <- function(x) dlnorm(x, meanlog[i], sdlog[i])^2
    integrate(d, 0, Inf)$value
  })
  
  return(c1 - c2)
}

# truncated-normal
qs.tnorm <- function(y, m, s, lb) {
  c1 <- 2*ftnorm(y, m, s, lb)
  c2 <- pnorm(lb, m, s, lower.tail = FALSE)^(-2) * (1 - pnorm(sqrt(2) * (lb - m)/s)) / (2 * s * sqrt(pi))

  return(c1 - c2)
}

################################################################################
### variable support

# generalized pareto distribution
qs.gpd <- function(y, location, scale, shape) {
  n <- max(lengths(list(y, location, scale, shape)))
  y <- rep(y, len = n)
  location <- rep(location, len = n)
  scale <- rep(scale, len = n)
  shape <- rep(shape, len = n)
  
  out <- numeric(n)
  ind <- abs(shape) > 1e-12
    
  upper <- ifelse(shape < 0, location - scale / shape, Inf)
  c1 <- 2*fgpd(y[ind], location[ind], scale[ind], shape[ind])
  c2 <- sapply(which(ind), function(i) {
    d <- function(x) fgpd(x, location[i], scale[i], shape[i])^2
    integrate(d, location[i], upper[i])$value
  })
  if (any(ind))
    warning("No closed form expression available for non-zero shape parameters - numerical integration used.")
  if (any(!ind & shape != 0))
    warning("Parameter 'shape' contains values close to zero. In those cases the QS is calculated assuming a value of 0.")
  out[ind] <- c1 - c2
  out[!ind] <- qs.exp(y - location, 1/scale)
  return(out)
}

# generalized extreme value distribution
qs.gev <- function(y, location, scale, shape) {
  n <- max(lengths(list(y, location, scale, shape)))
  y <- rep(y, len = n)
  location <- rep(location, len = n)
  scale <- rep(scale, len = n)
  shape <- rep(shape, len = n)
  
  ind1 <- abs(shape) > 1e-12
  ind2 <- shape > 0
  bound <- location - scale / shape
  lower <- ifelse(ind2, bound, -Inf)
  upper <- ifelse(ind2, Inf, bound)
  c1 <- 2*fgev(y, location, scale, shape)
  c2 <- numeric(n)
  
  c2[!ind1] <- 0.25/scale[!ind1]
  c2[ind1] <- sapply(which(ind1), function(i) {
    d <- function(x) fgev(x, location[i], scale[i], shape[i])^2
    integrate(d, lower[i], upper[i])$value
  })
  if (any(ind1))
    warning("No closed form expression available for non-zero shape parameters - numerical integration used.")
  if (any(!ind1 & shape != 0))
    warning("Parameter 'shape' contains values close to zero. In those cases the QS is calculated assuming a value of 0.")
  return(c1-c2)
}
