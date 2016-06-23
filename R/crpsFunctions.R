################################################################################
### general / non-parametric

# Numerical integration
# Note: y can be a vector, all other inputs are scalars
crps.int <- function(y, pxxx, lower, upper, rel_tol = 1e-6){
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

# (weighted) empirical distribution
crps.edf <- function(dat, y, w = NULL){
  n <- length(dat)
  # Set uniform weights unless specified otherwise
  if (is.null(w)){
    x <- sort(dat, decreasing = FALSE)
    out <- sapply(y, function(s) 2 / n^2 * sum((n * (s < x) - 1:n + 0.5) * (x - s)))
  } else {
    ord <- order(dat)
    x <- dat[ord]
    w <- w[ord]
    p <- c(0, cumsum(w))
    out <- sapply(y, function(s) 2 * sum((w * (s < x) - 0.5 * (p[2:(n+1)]^2 - p[1:n]^2)) * (x - s)))
  }
  return(out)
}

# kernel density estimation
crps.kdens = function(dat, y, bw = NULL){
  n <- length(dat)
  if (is.null(bw)) {
    s <- matrix(bw.nrd(dat), nrow = 1, ncol = n)
  }
  else {
    s <- matrix(bw, nrow = 1, ncol = n)
  }
  m <- matrix(dat, nrow = 1, ncol = n)
  w <- matrix(1/n, nrow = 1, ncol = n)
  return(crps.mixnorm(y = y, m = m, s = s, w = w))
}

################################################################################
### discrete / infinite support

# poisson
crps.pois <- function(y, lambda) {
  c1 <- (y - lambda) * (2*ppois(y, lambda) - 1)
  c2 <- 2*dpois(floor(y), lambda) - exp(-2*lambda) * (besselI(2*lambda, 0) + besselI(2*lambda, 1))
  return(c1 + lambda*c2)
}

# negative binomial
crps.nbinom <- function(y, size, prob) {
  if (!requireNamespace("hypergeo", quietly = TRUE)) {
    stop(paste(
      "Calculations require an implementation of the gaussian hypergeometric function.",
      "Please install the following package: hypergeo (>= 1.0)",
      sep = "\n"))
  }
  c1 <- y * (2 * pnbinom(y, size, prob) - 1)
  c2 <- (1 - prob) / prob ^ 2
  c3 <- prob * (2 * pnbinom(y - 1, size + 1, prob) - 1) + Re(hypergeo::hypergeo(size + 1, 0.5, 2,-4 * c2))
  return(c1 - size * c2 * c3)
}

################################################################################
### bounded interval

# uniform
crps.unif <- function(y, min, max) {
  c1 <- (y - min) * (2*punif(y, min, max) - 1)
  c2 <- (max - min) * (1/3 - punif(y, min, max)^2)
  return(c1 + c2)
}

# beta
crps.beta <- function(y, shape1, shape2) {
  c1 <- y * (2*pbeta(y, shape1, shape2) - 1)
  c2 <- shape1/(shape1+shape2)
  c3 <- 1 - 2*pbeta(y, shape1 + 1, shape2)
  c4 <- 2/shape1 * beta(2*shape1, 2*shape2) / beta(shape1, shape2)^2
  ind <- !is.finite(c4)
  if (any(ind)) {
    c4[ind] <- sqrt(shape2 / (pi * shape1 * (shape1 + shape2)))[ind]  # stirling's approximation
  }
  return(c1 + c2*(c3 - c4))
}

################################################################################
### real line

# laplace
crps.lapl <- function(y, location, scale) {
  z <- (y - location)/scale
  p <- 0.5 + 0.5 * sign(z) * pexp(abs(z))
  minp <- pmin(p, 1-p)
  c1 <- z*(2*p - 1) - 2*minp*(log(2*minp) - 1) - 0.75
  return(scale*c1)  
}

# logistic
crps.logis <- function(y, location, scale) {
  z <- (y - location)/scale
  p <- plogis(z)
  c1 <- z*(2*p - 1) - 1 - 2*(p*log(p) + (1-p)*log(1-p))
  return(scale*c1)
}

# normal
crps.norm <- function(y, location, scale,
                      lower = -Inf, upper = Inf,
                      lmass = 0, umass = 0) {
  
  ### standard formula
  
  ind1 <- any(is.finite(lower))
  ind2 <- any(is.finite(upper))
  
  if (!ind1 & !ind2) {
    z <- (y - location) / scale
    out_y <- z * (2 * pnorm(z) - 1) + 2 * dnorm(z) - 1 / sqrt(pi)
    return(scale * z)
  }
  
  ### dealing with truncation/censoring
  
  zb <- y
  if (ind1) {
    zb <- pmax(lower, zb)
    lb <- (lower - location) / scale
    
    if (is.character(lmass)) {
      Plb <- numeric(length(lmass))
      Plb[lmass == "cens"] <- pnorm(lb)
    } else {
      Plb <- lmass
    }
  }
  if (ind2) {
    zb <- pmin(upper, zb)
    ub <- (upper - location) / scale
    
    if (is.character(umass)) {
      Pub <- numeric(length(lmass))
      Pub[umass == "cens"] <- pnorm(ub, lower.tail = FALSE)
    } else {
      Pub <- umass
    }
  }
  res <- abs(y - zb)
  zb <- (zb - location) / scale
  
  if (ind1 & ind2) {
    a <- (1 - Plb - Pub) / (pnorm(ub) - pnorm(lb))
    out_l <- -lb * Plb^2 - 2 * a * dnorm(lb) * Plb + a^2 / sqrt(pi) * pnorm(lb * sqrt(2))
    out_u <- ub * Pub^2 - 2 * a * dnorm(ub) * Pub + a^2 / sqrt(pi) * pnorm(ub * sqrt(2), lower.tail = FALSE)
    out_y <- zb * (2 * (a * (pnorm(zb) - pnorm(lb)) + Plb) - 1) + 2 * a * dnorm(zb) - a^2 / sqrt(pi)
  } else if (ind1 & !ind2) {
    a <- (1 - Plb) / (1 - pnorm(lb))
    out_l <- -lb * Plb^2 - 2 * a * dnorm(lb) * Plb + a^2 / sqrt(pi) * pnorm(lb * sqrt(2))
    out_u <- 0
    out_y <- zb * (2 * (1 - a * pnorm(zb, lower.tail = FALSE)) - 1) + 2 * a * dnorm(zb) - a^2 / sqrt(pi)
  } else if (!ind1 & ind2) {
    a <- (1 - Pub) / pnorm(ub)
    out_l <- 0
    out_u <- ub * Pub^2 - 2 * a * dnorm(ub) * Pub + a^2 / sqrt(pi) * pnorm(ub * sqrt(2), lower.tail = FALSE)
    out_y <- zb * (2 * a * pnorm(zb) - 1) + 2 * a * dnorm(zb) - a^2 / sqrt(pi)
  }
  
  return(res + scale * (out_y + out_l + out_u))
}

# t
crps.t <- function(y, df, location, scale) {
  if (any(!df > 1))
    stop("Parameter 'df' contains values not greater than 1. The CRPS is not defined.")
  
  z <- (y - location) / scale
  ind <- df == Inf
  if (any(ind)) {
    if (length(z) < length(df)) {
      z <- rep(z, len = length(df))
    }
    out <- numeric(length(z))
    out[ind] <- crps.norm(z[ind], 0, 1)
    out[!ind] <- crps.t(z[!ind], df[!ind], 0, 1)
  } else {
    c1 <- z * (2 * pt(z, df) - 1)
    c2 <- dt(z, df) * (1 + z^2 / df)
    c3 <- beta(0.5, df - 0.5) / sqrt(df) / beta(0.5, 0.5 * df)^2
    out <- c1 + 2 * df / (df - 1) * (c2 - c3)
  }
  return(scale * out)
}

# mixture of normals (numerical integration)
crps.mixnorm.int <- function(y, m, s, w, rel_tol){
  Fmix <- function(z){
    sapply(z, function(r) sum(w*pnorm((r-m)/s)))
  }
  crps.int(y, Fmix, -Inf, Inf, rel_tol)
}

# mixture of normals
crps.mixnorm = function(y, m, s, w, exact = TRUE, rel_tol = 1e-6){
  if (exact == TRUE){
    out <- sapply(seq_along(y), function(i) crpsmixnC(w[i, ], m[i, ], s[i, ], y[i]))
  } else {
    out <- sapply(seq_along(y), function(i) crps.mixnorm.int(y[i], m[i, ], s[i, ], w[i, ], rel_tol))
  }
  return(out)
}

crps.2pexp <- function(y, location, scale1, scale2) {
  y1 <- pmin(y, location)
  y2 <- pmax(y, location)
  s <- scale1 + scale2
  a1 <- scale1 / s
  a2 <- scale2 / s
  b2 <- a1 - a2
  
  crps.cexp(-y1, -location, scale1, a1) +
    crps.cexp(y2, location, scale2, a2)
}

crps.2pnorm <- function(y, location, scale1, scale2) {
  y1 <- pmin(y, location)
  y2 <- pmax(y, location)
  s <- scale1 + scale2
  a1 <- scale1 / s
  a2 <- scale2 / s
  b2 <- a1 - a2
  
  crps.norm(y1, location, scale1, upper = location, umass = a2) +
    crps.norm(y2, location, scale2, lower = location, lmass = a1)
}

################################################################################
### non-negative

# exponential
crps.exp <- function(y, rate) {
  c1 <- 1/(2*rate) * (1 - 4*pexp(y, rate))
  return(abs(y) + c1)
}

# gamma
crps.gamma <- function(y, shape, scale) {
  c1 <- y*(2*pgamma(y, shape, scale=scale) - 1)
  c2 <- shape*(2*pgamma(y, shape+1, scale=scale) - 1) + 1/beta(.5, shape)
  return(c1 - scale*c2)
}

# log-laplace
crps.llapl <- function(y, locationlog, scalelog) {
  if (any(!scalelog < 1)) stop("Parameter 'scalelog' contains values not in (0, 1). The CRPS is not defined.")
  y1 <- log(pmax(y, 0))
  z <- (y1 - locationlog) / scalelog
  p <- 0.5 + 0.5 * sign(z) * pexp(abs(z))
  c1 <- y*(2*p - 1)
  c2 <- ifelse (z < 0,
          (1 - (2*p)^(1 + scalelog)) / (1 + scalelog),
          - (1 - (2*(1-p))^(1 - scalelog)) / (1 - scalelog)
  )
  c3 <- scalelog / (4 - scalelog^2) + c2
  return(c1 + exp(locationlog)*c3)
}

# log-logistic
crps.llogis <- function(y, locationlog, scalelog) {
  if (any(!scalelog < 1)) stop("Parameter 'scalelog' contains values not in (0, 1). The CRPS is not defined.")
  y1 <- log(pmax(y, 0))
  p <- plogis(y1, locationlog, scalelog)
  c1 <- y*(2*p - 1)
  c2 <- 2*exp(locationlog)*beta(1 + scalelog, 1 - scalelog)
  c3 <- (1 - scalelog)/2 - pbeta(p, 1 + scalelog, 1 - scalelog)
  return(c1 + c2*c3)
}

# log-normal
crps.lnorm <- function(y, meanlog, sdlog) {
  c1 <- y*(2*plnorm(y, meanlog, sdlog) - 1)
  c2 <- 2*exp(meanlog + 0.5*sdlog^2)
  c3 <- plnorm(y, meanlog + sdlog^2, sdlog) + pnorm(sdlog/sqrt(2)) - 1
  return(c1 - c2*c3)
}

# censored-exponential
crps.cexp <- function(y, location, scale, mass) {
  z <- (y - location)/scale
  c1 <- abs(z) - 2 * mass * pexp(z) + 0.5 * mass^2
  return(scale * c1)
}

################################################################################
### variable support

# generalized pareto distribution
crps.gpd <- function(y, location, scale, shape) {
  if (any(!shape < 1)) stop("Parameter 'shape' contains values not smaller than 1. The CRPS is not defined.")
  
  z <- (y - location)/scale
  ind <- abs(shape) < 1e-12
  if (any(ind)) {
    if (any(ind & shape != 0))
      warning("Parameter 'shape' contains values close to zero. In those cases the CRPS is calculated assuming a value of 0.")
    
    if (length(z) < length(shape)) {
      z <- rep(z, len = length(shape))
    }
    
    out <- numeric(length(z))
    out[ind] <- crps.exp(z[ind], 1)
    out[!ind] <- crps.gpd(z[!ind], 0, 1, shape[!ind])
  } else {
    p <- 1 - (1 + shape * z) ^ (-1 / shape)
    p[p < 0] <- 0
    p[p > 1] <- 1
    c1 <- (z + 1 / shape) * (2 * p - 1)
    c2 <- 2 / shape / (shape - 1) * (1 / (shape - 2) + (1 - p) ^ (1 - shape))
    out <- c1 - c2
  }
  return(scale * out)
}

# generalized extreme value distribution
crps.gev <- function(y, location, scale, shape) {
  if (any(!shape < 1)) stop("Parameter 'shape' contains values not smaller than 1. The CRPS is not defined.")
  
  z <- (y - location) / scale
  ind <- abs(shape) < 1e-12
  if (any(ind)) {
    if (any(ind & shape != 0))
      warning("Parameter 'shape' contains values close to zero. In those cases the CRPS is calculated assuming a value of 0.")
    
    if (length(z) < length(shape)) {
      z <- rep(z, len = length(shape))
    }
    
    out <- numeric(length(z))
    if (requireNamespace("gsl", quietly = TRUE)) {
      out[ind] <- -z[ind] - 2 * gsl::expint_Ei(-exp(-z[ind])) - digamma(1) - log(2)
    } else {
      warning(paste("The exponential integral is approximated using the 'integrate' function.",
                    "Consider installing the 'gsl' package to leverage a more accurate implementation.",
                    sep = "\n"))
      expint_Ei <- sapply(-exp(-z[ind]), function(upper) {
        integrate(function(x) exp(x)/x, -Inf, upper)$value
      })
      out[ind] <- -z[ind] - 2 * expint_Ei - digamma(1) - log(2)
    }
    out[!ind] <- crps.gev(z[!ind], 0, 1, shape[!ind])
  } else {
    x <- 1 + shape * z
    x[x < 0] <- 0
    p <- exp(-x^(-1/shape))
    out <- (-z - 1/shape)*(1 - 2*p) - 1/shape*gamma(1-shape)*(2^shape - 2*pgamma(-log(p), 1-shape))
  }
  return(scale * out)
}
