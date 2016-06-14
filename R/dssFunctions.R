################################################################################
### discrete / infinite support

# poisson
dss.pois <- function(y, lambda) {
  m <- lambda
  s <- sqrt(lambda)
  -dnorm(y, m, s, log=TRUE)
}

# negative binomial
dss.nbinom <- function(y, size, prob) { # hypergeo dependence
  m <- size * (1-prob)/prob
  s <- size * (1-prob)/prob^2
  -dnorm(y, m, s, log=TRUE)
}

################################################################################
### bounded interval

# uniform
dss.unif <- function(y, min, max) {
  m <- 0.5*(min + max)
  s <- 0.5*sqrt(1/3)*(max-min)
  -dnorm(y, m, s, log=TRUE)
}

# beta
dss.beta <- function(y, shape1, shape2) {
  m <- shape1 / (shape1 + shape2)
  s <- sqrt(shape1*shape2/(shape1 + shape2 + 1)) / (shape1 + shape2)
  -dnorm(y, m, s, log=TRUE)
}

################################################################################
### real line

# laplace
dss.lapl <- function(y, location, scale) {
  m <- location
  s <- sqrt(2) * scale
  -dnorm(y, m, s, log=TRUE)
}

# logistic
dss.logis <- function(y, location, scale) {
  m <- location
  s <- scale*pi/sqrt(3)
  -dnorm(y, m, s, log=TRUE)
}

# normal
dss.norm <- function(y, mean, sd) -dnorm(y, mean, sd, log=TRUE)

dss.mixn <- function(y, m, s, w) {
  mean <- rowSums(m * w)
  sd <- sqrt(rowSums((s^2 + (m - mean)^2) * w))
  -dnorm(y, mean, sd, log=TRUE)
}

# two-piece-normal
dss.2pnorm <- function(y, m, s1, s2) {
  mean <- m + sqrt(2/pi)*(s2 - s1)
  s <- sqrt((1-2/pi) * (s2 - s1)^2 + s1 * s2)
  -dnorm(y, mean, s, log=TRUE)
}

# t
dss.t <- function(y, df, location, scale) {
  if (any(!df>2)) stop("Parameter 'df' contains values not greater than 2. The DSS is not defined.")
  m <- location
  s <- scale * ifelse(is.finite(df), df / (df-2), 1)
  -dnorm(y, m, s, log=TRUE)
}


################################################################################
### non-negative

# exponential
dss.exp <- function(y, rate) {
  m <- 1/rate
  s <- 1/rate
  -dnorm(y, m, s, log=TRUE)
}

# gamma
dss.gamma <- function(y, shape, scale) {
  m <- shape * scale
  s <- sqrt(shape) * scale
  -dnorm(y, m, s, log=TRUE)
}

# log-laplace
dss.llapl <- function(y, locationlog, scalelog) {
  if (any(!scalelog < .5)) stop("Parameter 'scalelog contains values not in (0, 1/2). The DSS is not defined.")
  m <- exp(locationlog) / (1 - scalelog^2)
  s <- exp(locationlog) * sqrt(1 / (1 - 4*scalelog^2) - 1/(1-scalelog^2)^2)
  -dnorm(y, m, s, log=TRUE)
}

# log-logistic
dss.llogis <- function(y, locationlog, scalelog) {
  if (any(!scalelog < .5)) stop("Parameter 'scalelog contains values not in (0, 1/2). The DSS is not defined.")
  b <- pi*scalelog
  m <- exp(locationlog) * b / sin(b)
  s <- exp(locationlog) * sqrt(2*b/sin(2*b) - b^2/sin(b)^2)
  -dnorm(y, m, s, log=TRUE)
}

# log-normal
dss.lnorm <- function(y, meanlog, sdlog) {
  m <- exp(meanlog + sdlog^2/2)
  s <- sqrt(expm1(sdlog^2)*exp(2*meanlog + sdlog^2))
  -dnorm(y, m, s, log=TRUE)
}

# truncated-normal
dss.tn <- function(y, m, s, lb) {
  z <- (lb-m)/s
  mean <- m + s*dnorm(z)/(1-pnorm(z))
  sd <- s*sqrt((1 - dnorm(z)/(1-pnorm(z))*(dnorm(z)/(1-pnorm(z))-z)))
  -dnorm(y, mean, sd, log=TRUE)
}


################################################################################
### variable support

# generalized pareto distribution
dss.gpd <- function(y, location, scale, shape) {
  if (any(!shape<.5)) stop("Parameter 'shape' contains values not smaller than 1/2. The DSS does not exist.")
  m <- location + scale / (1 - shape)
  s <- scale / (1 - shape) / sqrt(1 - 2*shape)
  -dnorm(y, m, s, log=TRUE)
}

# generalized extreme value distribution
dss.gev <- function(y, location, scale, shape) {
  if (any(!shape < .5)) stop("Parameter 'shape' contains values not smaller than 1/2. The DSS does not exist.")
  ind <- abs(shape) < 1e-12
  m <- numeric(length(y))
  m[ind] <- location + scale * 0.57721566490153286
  m[!ind] <- location + scale / shape * (gamma(1-shape) - 1)
  s <- numeric(length(y))
  s[ind] <- scale * pi / sqrt(6)
  s[!ind] <- scale / shape * sqrt(gamma(1 - 2*shape) - gamma(1 - shape)^2)
  -dnorm(y, m, s, log=TRUE)
}