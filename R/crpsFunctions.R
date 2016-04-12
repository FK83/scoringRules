################################################################################
### general / non-parametric

# Numerical integration
crps.int <- function(y, pxxx, lower, upper){
  ind <- (y > upper) - (y < lower)
  out <- numeric(length(y))
  if (any(ind == -1)) {
    out[ind == -1] <- sapply(which(ind == -1), function(i) {
      F2 <- function(x) (1-pxxx[i](x))^2
      s1 <- lower[i] - y[i]
      s2 <- integrate(F2, lower[i], upper[i])$value
      s1 + s2
    })
  } else if (any(ind == 0)) {
    out[ind == 0] <- sapply(which(ind == 0), function(i) {
      F1 <- function(x) pxxx[i](x)^2
      F2 <- function(x) (1-pxxx[i](x))^2
      s1 <- integrate(F1, lower[i], y[i])$value
      s2 <- integrate(F2, y[i], upper[i])$value
      s1 + s2
    })
  } else if (any(ind == 1)) {
    out[ind == 1] <- sapply(which(ind == 1), function(i) {
      F1 <- function(x) pxxx[i](x)^2
      s1 <- integrate(F1, lower[i], upper[i])$value
      s2 <- y[i] - upper[i]
      s1 + s2
    })
  }
  return(-out)
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
  return(-out)
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
  return(-(c1 + lambda*c2))
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
  return(-(c1 - size * c2 * c3))
}

################################################################################
### bounded interval

# uniform
crps.unif <- function(y, min, max) {
  c1 <- (y - min) * (2*punif(y, min, max) - 1)
  c2 <- (max - min) * (1/3 - punif(y, min, max)^2)
  return(-(c1 + c2))
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
  return(-(c1 + c2*(c3 - c4)))
}

################################################################################
### real line

# laplace
crps.lapl <- function(y, location, scale) {
  z <- (y - location)/scale
  p <- 0.5 + 0.5 * sign(z) * pexp(abs(z))
  minp <- pmin(p, 1-p)
  c1 <- z*(2*p - 1) - 2*minp*(log(2*minp) - 1) - 0.75
  return(-scale*c1)  
}

# logistic
crps.logis <- function(y, location, scale) {
  z <- (y - location)/scale
  p <- plogis(z)
  c1 <- z*(2*p - 1) - 1 - 2*(p*log(p) + (1-p)*log(1-p))
  return(-scale*c1)
}

# normal
crps.norm <- function(y, mean, sd) {
  z <- (y - mean)/sd
  c1 <- z*(2*pnorm(z) - 1) + 2*dnorm(z) - 1/sqrt(pi)
  return(-sd*c1)
}

# mixture of normals
#crps.mixn = function (m, s, y, w = NULL, exact = TRUE, rel.tol = 1e-6){
#  n <- length(m)
#  if (is.null(w)) 
#    w <- rep(1/n, n)
#  if (exact == TRUE){
#    return(crpsmixnC(w, m, s, y))
#  } else {
#    Fmix = function(z){
#      sapply(z, function(r) sum(w*pnorm((r-m)/s)))
#    }
#    return(crps.int(Fmix, y, rel.tol = 1e-6))
# }
#}
crps.mixnorm <- function(y, m, s, w) {
  out <- sapply(seq_along(y), function(i) crpsmixnC(w[i, ], m[i, ], s[i, ], y[i]))
  return(out)
}

# two-piece-normal
crps.2pnorm <- function(y, m, s1, s2) {
  n <- max(lengths(list(y, m, s1, s2)))
  y <- rep(y, len = n)
  m <- rep(m, len = n)
  
  aux1 <- function(y, m, s1, s2) {
    a1 <- (y-m)/s1
    4*(s1^2)/(s1 + s2) * (a1*pnorm(a1) + dnorm(a1))
  }
  aux2 <- function(m, s1, s2) {
    2/(sqrt(pi)*(s1+s2)^2) * (sqrt(2)*s2*(s2^2-s1^2) - (s1^3+s2^3))    
  }
  sc <- ifelse(y <= m, aux1(y, m, s1, s2) - (y-m) + aux2(m, s1, s2), aux1(y, m, s2, s1) + (y-m)*((s1-s2)^2-4*s2^2)/(s1+s2)^2 + aux2(m, s2, s1))
  return(-sc)
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
    if (length(scale) < length(z)) {
      scale <- rep(scale, len = length(z))
    }
    out <- numeric(length(z))
    out[ind] <- scale[ind] * crps.norm(z[ind], 0, 1)
    out[!ind] <- scale[!ind] * crps.t(z[!ind], df[!ind], 0, 1)
  } else {
    c1 <- z * (2 * pt(z, df) - 1)
    c2 <- 2 * dt(z, df) * (df + z ^ 2) / (df - 1)
    c3 <- 2 * sqrt(df) / (df - 1) * beta(0.5, df - 0.5) / beta(0.5, 0.5 * df) ^ 2
    out <- -scale * (c1 + c2 - c3)
  }
  return(out)
}

################################################################################
### non-negative

# exponential
crps.exp <- function(y, rate) {
  c1 <- 1/(2*rate) * (1 - 4*pexp(y, rate))
  return(-(abs(y) + c1))
}

# gamma
crps.gamma <- function(y, shape, scale) {
  c1 <- y*(2*pgamma(y, shape, scale=scale) - 1)
  c2 <- shape*(2*pgamma(y, shape+1, scale=scale) - 1) - 1/beta(.5, shape)
  return(-(c1 - scale*c2))
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
  return(-(c1 + exp(locationlog)*c3))
}

# log-logistic
crps.llogis <- function(y, locationlog, scalelog) {
  if (any(!scalelog < 1)) stop("Parameter 'scalelog' contains values not in (0, 1). The CRPS is not defined.")
  y1 <- log(pmax(y, 0))
  p <- plogis(y1, locationlog, scalelog)
  c1 <- y*(2*p - 1)
  c2 <- 2*exp(locationlog)*beta(1 + scalelog, 1 - scalelog)
  c3 <- (1 - scalelog)/2 - pbeta(p, 1 + scalelog, 1 - scalelog)
  return(-(c1 + c2*c3))
}

# log-normal
crps.lnorm <- function(y, meanlog, sdlog) {
  c1 <- y*(2*plnorm(y, meanlog, sdlog) - 1)
  c2 <- 2*exp(meanlog + 0.5*sdlog^2)
  c3 <- plnorm(y, meanlog + sdlog^2, sdlog) + pnorm(sdlog/sqrt(2)) - 1
  return(-(c1 - c2*c3))
}

# truncated-normal
crps.tnorm <- function(y, m, s, lb) {
  z <- (y-m)/s
  x.lb <- (lb-m)/s; Phi.lb <- pnorm(-x.lb)
  return(-(s/Phi.lb^2*(z*Phi.lb*(2*(-pnorm(-z))+Phi.lb) + 2*dnorm(z)*Phi.lb - (pnorm(-sqrt(2)*x.lb))/sqrt(pi))))
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
    if (length(scale) < length(z)) {
      scale <- rep(scale, len = length(z))
    }
    
    out <- numeric(length(z))
    out[ind] <- scale[ind] * crps.exp(z[ind], 1)
    out[!ind] <- scale[!ind] * crps.gpd(z[!ind], 0, 1, shape[!ind])
  } else {
    p <- 1 - (1 + shape * z) ^ (-1 / shape)
    p[p < 0] <- 0
    p[p > 1] <- 1
    c1 <- (z + 1 / shape) * (2 * p - 1)
    c2 <- 2 / shape / (shape - 1) * (1 / (shape - 2) + (1 - p) ^ (1 - shape))
    out <- - scale * (c1 - c2)
  }
  return(out)
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
    if (length(scale) < length(z)) {
      scale <- rep(scale, len = length(scale))
    }
    
    out <- numeric(length(z))
    if (requireNamespace("gsl", quietly = TRUE)) {
      out[ind] <- -scale[ind] * (-z[ind] + (-digamma(1) - log(2)) - 2 * gsl::expint_Ei(-exp(-z[ind])))
    } else {
      warning(paste("Parameter 'shape' contains values equal or close to zero.",
                    "In those cases the CRPS is calculated using numerical integration and assuming a value of 0.",
                    "Consider installing the 'gsl' package to leverage an implementation of the exponential integral.",
                    sep = "\n"))
      if (length(ind) < length(z)) {
        ind <- rep(ind, len = length(z))
      }
      out[ind] <- sapply(which(ind), function(i) {
        F1 <- function(x) exp(-exp(-x))^2
        F2 <- function(x) (expm1(-exp(-x)))^2
        - scale[i] * (integrate(F1, -Inf, z[i])$value + integrate(F2, z[i], Inf)$value)
      })
    }
    out[!ind] <- scale[!ind] * crps.gev(z[!ind], 0, 1, shape[!ind])
  } else {
    p <- exp(-pmax(0, 1 + shape * z)^(-1/shape))
    out <- -scale * ((-z - 1/shape)*(1 - 2*p) - 1/shape*gamma(1-shape)*(2^shape - 2*pgamma(-log(p), 1-shape)))
  }
  return(out)
}
