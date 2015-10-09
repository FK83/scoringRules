################################################################################

list_qsFunctions <- list()

################################################################################
### bounded interval

# uniform
qs.unif <- function(y, min, max) 2*dunif(y, min, max) - 1/(max-min)
list_qsFunctions$'uniform' <- "qs.unif"

# beta
qs.beta <- function(y, shape1, shape2) {
  c1 <- 2 * dbeta(y, shape1, shape2)
  c2 <- beta(2*shape1 - 1, 2*shape2 - 1) / beta(shape1, shape2)^2
  return(c1 - c2)
}
list_qsFunctions$'beta' <- "qs.beta"


################################################################################
### real line

# laplace
qs.lapl <- function(y, location, scale) {
  c1 <- 2*flapl(y, location, scale)
  c2 <- 1/(4*b)
  return(c1 - c2)
}
list_qsFunctions$'laplace' <- "qs.lapl"

# logistic
qs.logis <- function(y, location, scale) {
  c1 <- 2*dlogis(y, location, scale)
  c2 <- sapply(seq_along(y), function(i) {
    d <- function(x) dlogis(x, location[i], scale[i])^2
    integrate(d, -Inf, Inf)$value
  })
  return(c1 - c2)
}
list_qsFunctions$'logistic' <- "qs.logis"

# normal
qs.norm <- function(y, mean, sd) {
  c1 <- 2*dnorm(y, mean, sd)
  c2 <- 1/(2*sd*sqrt(pi))
  return(c1 - c2)
}
list_qsFunctions$'normal' <- "qs.norm"

# mixture of normals
qs.mixn <- function(y, m, s, w) sapply(seq_along(y), function(i) qsmixnC(w[i, ], m[i, ], s[i, ], y[i]))
list_qsFunctions$'normal-mixture' <- "qs.mixn"

# two-piece-normal
qs.2pnorm <- function(y, m, s1, s2) {
  c1 <- 2*f2pnorm(y, m, s1, s2)
  c2 <- (s1 + s2) / (sqrt(pi) * (s1 + s2)^2)
  return(c1 - c2)
}
list_qsFunctions$'two-piece-normal' <- "qs.2pnorm"

# t
qs.t <- function(y, df, location, scale) {
  c1 <- 2*ft(y, df, location, scale)
#  c2 <- 1/df^(1.5) / sqrt(pi) / scale * gamma(df + 0.5) / gamma(df) * gamma((df + 1)/2)^2 / gamma(df/2)^2
#  c3 <- 2^(3-4*df) * pi / scale / df^(1.5) * gamma(2*df) / gamma(df/2)^4
  c2 <- 2^(3-4*df) * pi / scale / df^(1.5) / beta(df, df) / beta(df/2, df/2)^2
  ind <- !is.finite(c2)
  c2[ind] <- 1/(2*scale*sqrt(pi))
  return(c1 - c2)
}
list_qsFunctions$'t' <- "qs.t"


################################################################################
### non-negative

# exponential
qs.exp <- function(y, rate) 2*dexp(y, rate) - 0.5*rate
list_qsFunctions$'exponential' <- "qs.exp"

# gamma
qs.gamma <- function(y, shape, scale) {
  c1 <- 2*dgamma(y, shape, scale = scale)
  c2 <- 2^(1-2*shape) / scale / (2*shape - 1) / beta(shape, shape)
  return(c1 - c2)
}
list_qsFunctions$'gamma' <- "qs.gamma"

# log-laplace
qs.llapl <- function(y, locationlog, scalelog) {
  c1 <- 2*fllapl(y, locationlog, scalelog)
  c2 <- sapply(seq_along(y), function(i) {
    d <- function(x) fllapl(x, locationlog[i], scalelog[i])^2
    integrate(d, 0, Inf)$value
  }) 
  return(c1 - c2)
}
list_qsFunctions$'log-laplace' <- "qs.llapl"

# log-logistic
qs.llogis <- function(y, locationlog, scalelog) {
  c1 <- 2*fllogis(y, locationlog, scalelog)
  c2 <- sapply(seq_along(y), function(i) {
    d <- function(x) fllogis(x, locationlog[i], scalelog[i])^2
    integrate(d, 0, Inf)$value
  })
  return(c1 - c2)
}
list_qsFunctions$'log-logistic' <- "qs.llogis"

# log-normal
qs.lnorm <- function(y, meanlog, sdlog) {
  c1 <- 2*dlnorm(y, meanlog, sdlog)
  c2 <- sapply(seq_along(y), function(i) {
    d <- function(x) dlnorm(x, meanlog[i], sdlog[i])^2
    integrate(d, 0, Inf)$value
  })
  return(c1 - c2)
}
list_qsFunctions$'log-normal' <- "qs.lnorm"

# truncated-normal
qs.tn <- function(y, m, s, lb) {
  c1 <- 2*ftn(y, m, s, lb)
  c2 <- sapply(seq_along(y), function(i) {
    d <- function(x) ftn(x, m[i], s[i], lb[i])^2
    integrate(d, lb[i], Inf)$value
  })
  return(c1 - c2)
}
list_qsFunctions$'truncated-normal' <- "qs.tn"


################################################################################
### variable support

# generalized pareto distribution
qs.gpd <- function(y, location, scale, shape) {
  ind <- abs(shape) > 1e-12
  out <- numeric(length(y))
    
  if (any(!ind)) {
    out[!ind] <- qs.exp(y - location, 1/scale)
  }
    
  upper <- ifelse(shape < 0, location - scale / shape, Inf)
  c1 <- 2*fgpd(y[ind], location[ind], scale[ind], shape[ind])
  c2 <- sapply(which(ind), function(i) {
    d <- function(x) fgpd(x, location[i], scale[i], shape[i])^2
    integrate(d, location[i], upper[i])$value
  })
  out[ind] <- c1 - c2
  return(out)
}
list_qsFunctions$'gpd' <- "qs.gpd"

# generalized extreme value distribution
qs.gev <- function(y, location, scale, shape) {
  ind1 <- abs(shape) > 1e-12
  ind2 <- shape > 0
  bound <- location - scale / shape
  lower <- ifelse(ind2, bound, -Inf)
  upper <- ifelse(ind2, Inf, bound)
  c1 <- 2*fgev(y, location, scale, shape)
  c2 <- numeric(length(y))
  
  c2[!ind1] <- 0.25/scale[!ind1]
  c2[ind1] <- sapply(which(ind1), function(i) {
    d <- function(x) fgev(x, location[i], scale[i], shape[i])^2
    integrate(d, lower[i], upper[i])$value
  })
  
  return(c1-c2)
}
list_qsFunctions$'gev' <- "qs.gev"
