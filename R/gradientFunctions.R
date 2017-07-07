

### gradient (location, scale) for censored normal distribution
gradcnorm <- function(y, location, scale,
                      lower = -Inf, upper = Inf) {
  
  z <- (pmin(pmax(y, lower), upper) - location) / scale
  lb <- (lower - location) / scale
  ub <- (upper - location) / scale
  
  term0 <- crps.norm(z, 0, 1, lb, ub, "cens", "cens")
  term1 <- 2 * pnorm(z) - 1
  term2 <- pnorm(lb)^2
  term3 <- pnorm(ub, lower.tail = FALSE)^2
  
  dmu <- -term1 + term2 - term3
  dsigma <- term0 - term1 * z +
    term2 * ifelse(is.finite(lb), lb, 0) -
    term3 * ifelse(is.finite(ub), ub, 0)
  
  return(cbind(dmu, dsigma))
}

### gradient (location, scale) for truncated normal distribution
gradtnorm <- function(y, location, scale,
                      lower = -Inf, upper = Inf) {
  
  z <- (pmin(pmax(y, lower), upper) - location) / scale
  lb <- (lower - location) / scale
  ub <- (upper - location) / scale
  
  renorm <- 1 / (pnorm(ub) - pnorm(lb))
  FF <- renorm * (pnorm(z) - pnorm(lb))
  
  term0 <- crps.norm(z, 0, 1, lb, ub, "trunc", "trunc")
  term1 <- 2 * FF  - 1
  term2 <- 2 * renorm * dnorm(lb) *
    (z * FF + renorm * (dnorm(z) - dnorm(lb)) - term0)
  term3 <- 2 * renorm * dnorm(ub) *
    (z * (1 - FF) + renorm * (dnorm(ub) - dnorm(z)) + term0)
  
  dmu <- -term1 + term2 + term3
  dsigma <- term0 - term1 * z + 
    term2 * ifelse(is.finite(lb), lb, 0) +
    term3 * ifelse(is.finite(ub), ub, 0)
  
  return(cbind(dmu, dsigma))
}

### gradient (location, scale) for censored t distribution
gradct <- function(y, df, location, scale,
                   lower = -Inf, upper = Inf) {
  
  z <- (pmin(pmax(y, lower), upper) - location) / scale
  lb <- (lower - location) / scale
  ub <- (upper - location) / scale
  
  term0 <- crps.t(z, df, 0, 1, lb, ub, "cens", "cens")
  term1 <- 2 * pt(z, df) - 1
  term2 <- pt(lb, df)^2
  term3 <- pt(ub, df, lower.tail = FALSE)^2
  
  dmu <- -term1 + term2 - term3
  dsigma <- term0 - term1 * z +
    term2 * ifelse(is.finite(lb), lb, 0) -
    term3 * ifelse(is.finite(ub), ub, 0)
  
  return(cbind(dmu, dsigma))
}

### gradient (location, scale) for truncated t distribution
gradtt <- function(y, df, location, scale,
                   lower = -Inf, upper = Inf) {
  
  z <- (pmin(pmax(y, lower), upper) - location) / scale
  lb <- (lower - location) / scale
  ub <- (upper - location) / scale
  
  renorm <- 1 / (pt(ub, df) - pt(lb, df))
  FF <- renorm * (pt(z, df) - pt(lb, df))
  G_renorm <- renorm * ifelse(is.finite(df), df / (df - 1), 1)
  Gz <- (1 + z^2/df) * dt(z, df)
  Glb <- ifelse(is.finite(lb), (1 + lb^2/df) * dt(lb, df), 0)
  Gub <- ifelse(is.finite(ub), (1 + ub^2/df) * dt(ub, df), 0)
  
  term0 <- crps.t(z, df, 0, 1, lb, ub, "trunc", "trunc")
  term1 <- 2 * FF  - 1
  term2 <- 2 * renorm * dt(lb, df) *
    (z * FF + G_renorm * (Gz - Glb) - term0)
  term3 <- 2 * renorm * dt(ub, df) *
    (z * (1 - FF) + G_renorm * (Gub - Gz) + term0)
  
  dmu <- -term1 + term2 + term3
  dsigma <- term0 - term1 * z + 
    term2 * ifelse(is.finite(lb), lb, 0) +
    term3 * ifelse(is.finite(ub), ub, 0)
  
  return(cbind(dmu, dsigma))
}

### gradient (location, scale) for censored logistic distribution
gradclogis <- function(y, location, scale,
                       lower = -Inf, upper = Inf) {
  
  z <- (pmin(pmax(y, lower), upper) - location) / scale
  lb <- (lower - location) / scale
  ub <- (upper - location) / scale
  
  term0 <- crps.logis(z, 0, 1, lb, ub, "cens", "cens")
  term1 <- 2 * plogis(z) - 1
  term2 <- plogis(lb)^2
  term3 <- plogis(ub, lower.tail = FALSE)^2
  
  dmu <- -term1 + term2 - term3
  dsigma <- term0 - term1 * z +
    term2 * ifelse(is.finite(lb), lb, 0) -
    term3 * ifelse(is.finite(ub), ub, 0)
  
  return(cbind(dmu, dsigma))
}

### gradient (location, scale) for truncated logistic distribution
gradtlogis <- function(y, location, scale,
                       lower = -Inf, upper = Inf) {
  
  z <- (pmin(pmax(y, lower), upper) - location) / scale
  lb <- (lower - location) / scale
  ub <- (upper - location) / scale
  
  renorm <- 1 / (plogis(ub) - plogis(lb))
  FF <- renorm * (plogis(z) - plogis(lb))
  Gz <- plogis(z, log.p = TRUE) - z * plogis(z, lower.tail = FALSE)
  Glb <- ifelse(
    is.finite(lb),
    plogis(lb, log.p = TRUE, lower.tail = FALSE) + lb * plogis(lb),
    0
  )
  Gub <- ifelse(
    is.finite(ub),
    plogis(ub, log.p = TRUE) - ub * plogis(ub, lower.tail = FALSE),
    0
  )
  
  term0 <- crps.logis(z, 0, 1, lb, ub, "trunc", "trunc")
  term1 <- 2 * FF  - 1
  term2 <- 2 * renorm * dlogis(lb) *
    (z * FF - renorm * (Gz - Glb) - term0)
  term3 <- 2 * renorm * dlogis(ub) *
    (z * (1 - FF) - renorm * (Gub - Gz) + term0)
  
  dmu <- -term1 + term2 + term3
  dsigma <- term0 - term1 * z + 
    term2 * ifelse(is.finite(lb), lb, 0) +
    term3 * ifelse(is.finite(ub), ub, 0)
  
  return(cbind(dmu, dsigma))
}
