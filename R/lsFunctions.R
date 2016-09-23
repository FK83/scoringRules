################################################################################
### discrete / infinite support

# poisson
ls.pois <- function(y, lambda) -dpois(y, lambda, log=TRUE)

# negative binomial
ls.nbinom <- function(y, size, prob) -dnbinom(y, size, prob, log=TRUE)

################################################################################
### bounded interval

# uniform
ls.unif <- function(y, min, max, lmass = 0, umass = 0) {
  if (any(lmass != 0 | umass != 0)) stop("Log score unavailable due to point mass.")
  -dunif(y, min, max, log=TRUE)
}

# beta
ls.beta <- function(y, shape1, shape2) -dbeta(y, shape1, shape2, log=TRUE)

################################################################################
### real line

# laplace
ls.lapl <- function(y, location, scale) -log(flapl(y, location, scale))

# logistic
ls.logis <- function(y, location, scale,
                     lower = -Inf, upper = Inf,
                     lmass = 0, umass = 0) {
  if (any(lmass != 0 | umass != 0)) stop("Log score unavailable due to point mass.")
  -flogis(y, location, scale, lower, upper, lmass, umass, log = TRUE)
}

# normal
ls.norm <- function(y, location, scale,
                    lower = -Inf, upper = Inf,
                    lmass = 0, umass = 0) {
  if (any(lmass != 0 | umass != 0)) stop("Log score unavailable due to point mass.")
  -fnorm(y, location, scale, lower, upper, lmass, umass, log = TRUE)
}

# t
ls.t <- function(y, df, location, scale,
                 lower = -Inf, upper = Inf,
                 lmass = 0, umass = 0) {
  if (any(lmass != 0 | umass != 0)) stop("Log score unavailable due to point mass.")
  -ft(y, df, location, scale, lower, upper, lmass, umass, log = TRUE)
}

# mixture of normals
ls.mixnorm <- function(y, m, s, w) sapply(seq_along(y), function(i) lsmixnC(w[i, ], m[i, ], s[i, ], y[i]))

# two-piece-exponential
ls.2pexp <- function(y, location, scale1, scale2) -log(f2pexp(y, location, scale1, scale2))

# two-piece-normal
ls.2pnorm <- function(y, location, scale1, scale2) -log(f2pnorm(y, location, scale1, scale2))

################################################################################
### non-negative

# gamma
ls.gamma <- function(y, shape, scale) -dgamma(y, shape, scale = scale, log=TRUE)

# log-laplace
ls.llapl <- function(y, locationlog, scalelog) -log(fllapl(y, locationlog, scalelog))

# log-logistic
ls.llogis <- function(y, locationlog, scalelog) -log(fllogis(y, locationlog, scalelog))

# log-normal
ls.lnorm <- function(y, meanlog, sdlog) -dlnorm(y, meanlog, sdlog, log=TRUE)

# censored-exponential
ls.cexp <- function(y, location, scale, mass) {
  stop("")
}

################################################################################
### variable support

# exponential
ls.exp <- function(y, location, scale, mass = 0) {
  if (any(mass != 0)) stop("Log score unavailable due to point mass.")
  -fexp(y, location, scale, mass, log = TRUE)
}

# generalized pareto distribution
ls.gpd <- function(y, location, scale, shape, mass = 0) {
  if (any(mass != 0)) stop("Log score unavailable due to point mass.")
  -fgpd(y, location, scale, shape, mass, log = TRUE)
} 

# generalized extreme value distribution
ls.gev <- function(y, location, scale, shape) -log(fgev(y, location, scale, shape))
