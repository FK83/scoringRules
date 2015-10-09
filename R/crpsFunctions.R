################################################################################

list_crpsFunctions <- list()

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
    x <- .Internal(sort(dat, FALSE))
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
  return(crps.mixn(y = y, m = m, s = s, w = w))
}

################################################################################
### discrete / infinite support

# poisson
crps.pois <- function(y, lambda) {
  c1 <- (y - lambda) * (2*ppois(y, lambda) - 1)
  c2 <- 2*dpois(floor(y), lambda) - exp(-2*lambda) * (besselI(2*lambda, 0) + besselI(2*lambda, 1))
  return(-(c1 + lambda*c2))
}
list_crpsFunctions$'poisson' <- "crps.pois"

# negative binomial
crps.nbinom <- function(y, size, prob) { # hypergeo dependence
  c1 <- y*(2*pnbinom(y, size, prob) - 1)
  c2 <- (1-prob)/prob^2
  c3 <- prob*(2*pnbinom(y-1, size+1, prob) - 1) + as.numeric(hypergeo(size+1, 0.5, 2, -4*c2))
  return(-(c1 - size*c2*c3))
}
list_crpsFunctions$'negative-binomial' <- "crps.nbinom"

################################################################################
### bounded interval

# uniform
crps.unif <- function(y, min, max) {
  c1 <- (y - min) * (2*punif(y, min, max) - 1)
  c2 <- (max - min) * (1/3 - punif(y, min, max)^2)
  return(-(c1 + c2))
}
list_crpsFunctions$'uniform' <- "crps.unif"

# beta
crps.beta <- function(y, shape1, shape2) {
  c1 <- y * (2*pbeta(y, shape1, shape2) - 1)
  c2 <- shape1/(shape1+shape2)
  c3 <- 1 - 2*pbeta(y, shape1 + 1, shape2)
  c4 <- 2/shape1 * beta(2*shape1, 2*shape2) / beta(shape1, shape2)^2
  ind <- !is.finite(c4)
  c4[ind] <- sqrt(shape2[ind] / (pi*shape1[ind]*(shape1[ind]+shape2[ind])))
  return(-(c1 + c2*(c3 - c4)))
}
list_crpsFunctions$'beta' <- "crps.beta"


################################################################################
### real line

# laplace
crps.lapl <- function(y, location, scale) {
  z <- (y - location)/scale
  p <- ifelse(z < 0, 0.5 * exp(z), 1 - 0.5 * exp(-z))
  minp <- pmin(p, 1-p)
  c1 <- z*(2*p - 1) - 2*minp*(log(2*minp) - 1) - 0.75
  return(-scale*c1)  
}
list_crpsFunctions$'laplace' <- "crps.lapl"

# logistic
crps.logis <- function(y, location, scale) {
  z <- (y - location)/scale
  p <- plogis(z)
  c1 <- z*(2*p - 1) - 1 - 2*(p*log(p) + (1-p)*log(1-p))
  return(-scale*c1)
}
list_crpsFunctions$'logistic' <- "crps.logis"

# normal
crps.norm <- function(y, mean, sd) {
  z <- (y - mean)/sd
  c1 <- z*(2*pnorm(z) - 1) + 2*dnorm(z) - 1/sqrt(pi)
  return(-sd*c1)
}
list_crpsFunctions$'normal' <- "crps.norm"

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
crps.mixn <- function(y, m, s, w) {
  out <- sapply(seq_along(y), function(i) crpsmixnC(w[i, ], m[i, ], s[i, ], y[i]))
  return(out)
}
list_crpsFunctions$'normal-mixture' <- "crps.mixn"

# two-piece-normal
crps.2pnorm <- function(y, m, s1, s2) {
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
list_crpsFunctions$'two-piece-normal' <- "crps.2pnorm"

# t
crps.t <- function(y, df, location, scale) {
  if (any(!df > 1)) stop("Parameter 'df' contains values not greater than 1. The CRPS is not defined.")
  ifelse (df == Inf,
          return(crps.norm(y, location, scale)),
          {
            z <- (y - location)/scale
            c1 <- z*(2*pt(z, df) - 1)
            c2 <- 2*dt(z, df) * (df + z^2)/(df - 1)
            c3 <- 2*sqrt(df)/(df - 1) * beta(0.5, df - 0.5) / beta(0.5, 0.5*df)^2
            return(-scale*(c1 + c2 - c3))
          }
  )
}
list_crpsFunctions$'t' <- "crps.t"


################################################################################
### non-negative

# exponential
crps.exp <- function(y, rate) {
  c1 <- 1/(2*rate) * (1 - 4*pexp(y, rate))
  return(-(abs(y) + c1))
}
list_crpsFunctions$'exponential' <- "crps.exp"

# gamma
crps.gamma <- function(y, shape, scale) {
  c1 <- y*(2*pgamma(y, shape, scale=scale) - 1)
  c2 <- shape*(2*pgamma(y, shape+1, scale=scale) - 1) - 1/beta(.5, shape)
  return(-(c1 - scale*c2))
}
list_crpsFunctions$'gamma' <- "crps.gamma"

# log-laplace
crps.llapl <- function(y, locationlog, scalelog) {
  if (any(!scalelog < 1)) stop("Parameter 'scalelog' contains values not in (0, 1). The CRPS is not defined.")
  y1 <- suppressWarnings(ifelse(y > 0, log(y), -Inf))
  z <- (y1 - locationlog) / scalelog
  p <- ifelse(z < 0, 0.5 * exp(z), 1 - 0.5 * exp(-z))
  scale <- exp(locationlog)
  c1 <- y*(2*p - 1)
  c2 <- ifelse (y < scale,
          (1 - (2*p)^(1 + scalelog)) / (1 + scalelog),
          - (1 - (2*(1-p))^(1 - scalelog)) / (1 - scalelog)
  )
  c3 <- scalelog / (4 - scalelog^2) + c2
  return(-(c1 + scale*c3))
}
list_crpsFunctions$'log-laplace' <- "crps.llapl"

# log-logistic
crps.llogis <- function(y, locationlog, scalelog) {
  if (any(!scalelog < 1)) stop("Parameter 'scalelog' contains values not in (0, 1). The CRPS is not defined.")
  y1 <- suppressWarnings(ifelse(y > 0, log(y), -Inf))
  p <- plogis(y1, locationlog, scalelog)
  c1 <- y*(2*p - 1)
  c2 <- 2*exp(locationlog)*beta(1 + scalelog, 1 - scalelog)
  c3 <- (1 - scalelog)/2 - pbeta(p, 1 + scalelog, 1 - scalelog)
  return(-(c1 + c2*c3))
}
list_crpsFunctions$'log-logistic' <- "crps.llogis"

# log-normal
crps.lnorm <- function(y, meanlog, sdlog) {
  c1 <- y*(2*plnorm(y, meanlog, sdlog) - 1)
  c2 <- 2*exp(meanlog + 0.5*sdlog^2)
  c3 <- plnorm(y, meanlog + sdlog^2, sdlog) + pnorm(sdlog/sqrt(2)) - 1
  return(-(c1 - c2*c3))
}
list_crpsFunctions$'log-normal' <- "crps.lnorm"

# truncated-normal
crps.tn <- function(y, m, s, lb) {
  z <- (y-m)/s
  x.lb <- (lb-m)/s; Phi.lb <- pnorm(-x.lb)
  return(-(s/Phi.lb^2*(z*Phi.lb*(2*(-pnorm(-z))+Phi.lb) + 2*dnorm(z)*Phi.lb - (pnorm(-sqrt(2)*x.lb))/sqrt(pi))))
}
list_crpsFunctions$'truncated-normal' <- "crps.tn"


################################################################################
### variable support

# generalized pareto distribution
crps.gpd <- function(y, location, scale, shape) {
  if (any(!shape < 1)) stop("Parameter 'shape' contains values not smaller than 1. The CRPS is not defined.")
  ind <- abs(shape) < 1e-12
  out <- numeric(length(y))
  
  if (any(ind)) {
    if (any(ind & shape != 0)) warning("Parameter 'shape' contains values close to zero. In those cases the CRPS is calculated assuming a value of 0.")
    out[ind] <- crps.exp(y[ind] - location[ind], 1/scale[ind])
    y <- y[!ind]
    location <- location[!ind]
    scale <- scale[!ind]
    shape <- shape[!ind]
  }
  
  z <- (y - location) / scale
  p <- pmin(1, ifelse(y < location, 0, 1-(1 + shape * z)^(-1/shape)))
  c1 <- (y - location + scale/shape) * (2*p1 - 1)
  c2 <- 2*scale/shape/(shape-1) * (1/(shape-2) + (1-p)^(1-shape))
  out[!ind] <- -(c1 - c2)
  
  return(out)
}
list_crpsFunctions$'gpd' <- "crps.gpd"

# generalized extreme value distribution
crps.gev <- function(y, location, scale, shape) {
  if (any(!shape < 1)) stop("Parameter 'shape' contains values not smaller than 1. The CRPS is not defined.")
  ind <- abs(shape) < 1e-12
  out <- numeric(length(y))
  
  if (any(ind)) {
    warning("Parameter 'shape' contains values equal or close to zero.\n  In those cases the CRPS is calculated using numerical integration and assuming a value of 0.")
    out[ind] <- sapply(which(ind), 1, function(i) {
      F1 <- function(x) exp(-exp(-(x-location[i])/scale[i]))^2
      F2 <- function(x) (expm1(-exp(-(x-location[i])/scale[i])))^2
      integrate(F1, -Inf, y[i])$value + integrate(F2, y[i], Inf)$value
    })
    y <- y[!ind]
    location <- location[!ind]
    scale <- scale[!ind]
    shape <- shape[!ind]
  }
  
  z <- (y - location) / scale
  p <- exp(-pmax(0, 1 + shape * z)^(-1/shape))
  out[!ind] <- (location - y - scale/shape)*(1 - 2*p) - scale/shape*gamma(1-shape)*(2^shape - 2*pgamma(-log(p), 1-shape))
  
  return(-out)
}
list_crpsFunctions$'gev' <- "crps.gev"

