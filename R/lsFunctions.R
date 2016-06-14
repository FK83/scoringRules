################################################################################
### bounded interval

# uniform
ls.unif <- function(y, min, max) -dunif(y, min, max, log=TRUE)

# beta
ls.beta <- function(y, shape1, shape2) -dbeta(y, shape1, shape2, log=TRUE)

################################################################################
### real line

# laplace
ls.lapl <- function(y, location, scale) -log(flapl(y, location, scale))

# logistic
ls.logis <- function(y, location, scale) -dlogis(y, location, scale, log=TRUE)

# normal
ls.norm <- function(y, mean, sd) -dnorm(y, mean, sd, log=TRUE)

# mixture of normals
ls.mixnorm <- function(y, m, s, w) -sapply(seq_along(y), function(i) lsmixnC(w[i, ], m[i, ], s[i, ], y[i]))

# two-piece-normal
ls.2pnorm <- function(y, m, s1, s2) -log(f2pnorm(y, m, s1, s2))

# t
ls.t <- function(y, df, location, scale) -log(ft(y, df, location, scale))

################################################################################
### non-negative

# exponential
ls.exp <- function(y, rate) -dexp(y, rate, log=TRUE)

# gamma
ls.gamma <- function(y, shape, scale) -dgamma(y, shape, scale = scale, log=TRUE)

# log-laplace
ls.llapl <- function(y, locationlog, scalelog) -log(fllapl(y, locationlog, scalelog))

# log-logistic
ls.llogis <- function(y, locationlog, scalelog) -log(fllogis(y, locationlog, scalelog))

# log-normal
ls.lnorm <- function(y, meanlog, sdlog) -dlnorm(y, meanlog, sdlog, log=TRUE)

# truncated-normal
ls.tnorm <- function(y, m, s, lower, upper) -log(ftnorm(y, m, s, lower, upper))

################################################################################
### variable support

# generalized pareto distribution
ls.gpd <- function(y, location, scale, shape) -log(fgpd(y, location, scale, shape))

# generalized extreme value distribution
ls.gev <- function(y, location, scale, shape) -log(fgev(y, location, scale, shape))
