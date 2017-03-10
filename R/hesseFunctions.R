### Hessian (location, scale) for censored normal distribution
hesscnorm <- function(y, location, scale,
                      lower = -Inf, upper = Inf) {
  
  z <- (pmin(pmax(y, lower), upper) - location) / scale
  lb <- (lower - location) / scale
  ub <- (upper - location) / scale
  
  term1 <- dnorm(z)
  term2 <- dnorm(lb) * pnorm(lb)
  term3 <- dnorm(ub) * pnorm(ub, lower.tail = FALSE)
  
  d2mu <- term1 - term2 - term3
  dmu.dsigma <- dsigma.dmu <- term1 * z -
    term2 * ifelse(is.finite(lb), lb, 0) -
    term3 * ifelse(is.finite(ub), ub, 0)
  d2sigma <- term1 * z^2 -
    term2 * ifelse(is.finite(lb), lb^2, 0) -
    term3 * ifelse(is.finite(ub), ub^2, 0)
  
  return(2/scale * cbind(d2mu, d2sigma, dmu.dsigma, dsigma.dmu))
}

### Hessian (location, scale) for truncated normal distribution
hesstnorm <- function(y, location, scale,
                      lower = -Inf, upper = Inf) {
  
  z <- (pmin(pmax(y, lower), upper) - location) / scale
  lb <- (lower - location) / scale
  ub <- (upper - location) / scale
  z_lb_ub <- cbind(z, lb, ub)
  
  CDF <- pnorm(z_lb_ub)
  denom <- CDF[, "ub"] - CDF[, "lb"]
  CDF <- CDF / denom
  PDF <- dnorm(z_lb_ub) / denom
  PDFp_over_PDF <- -z_lb_ub
  G <- -PDF
  CRPS <- crps.norm(z, 0, 1, lb, ub, "trunc", "trunc")
  
  Hessian <- calcHess_trunc(z, scale, lb, ub,
                            CDF, PDF, PDFp_over_PDF, G, CRPS)
  
  return(Hessian)
}

### Hessian (location, scale) for censored t distribution
hessct <- function(y, df, location, scale,
                      lower = -Inf, upper = Inf) {
  
  z <- (pmin(pmax(y, lower), upper) - location) / scale
  lb <- (lower - location) / scale
  ub <- (upper - location) / scale
  
  term1 <- dt(z, df)
  term2 <- dt(lb, df) * pt(lb, df)
  term3 <- dt(ub, df) * pt(ub, df, lower.tail = FALSE)
  
  d2mu <- term1 - term2 - term3
  dmu.dsigma <- dsigma.dmu <- term1 * z -
    term2 * ifelse(is.finite(lb), lb, 0) -
    term3 * ifelse(is.finite(ub), ub, 0)
  d2sigma <- term1 * z^2 -
    term2 * ifelse(is.finite(lb), lb^2, 0) -
    term3 * ifelse(is.finite(ub), ub^2, 0)
  
  return(2/scale * cbind(d2mu, d2sigma, dmu.dsigma, dsigma.dmu))
}

### gradient (location, scale) for truncated t distribution
hesstt <- function(y, df, location, scale,
                   lower = -Inf, upper = Inf) {
  
  z <- (pmin(pmax(y, lower), upper) - location) / scale
  lb <- (lower - location) / scale
  ub <- (upper - location) / scale
  z_lb_ub <- cbind(z, lb, ub)
  
  CDF <- pt(z_lb_ub, df)
  denom <- CDF[, "ub"] - CDF[, "lb"]
  CDF <- CDF / denom
  PDF <- dt(z_lb_ub, df) / denom
  PDFp_over_PDF <- -(df + 1) / df / (1 + z_lb_ub^2/df) *
    ifelse(is.finite(z_lb_ub), z_lb_ub, 0)
  G <- -df / (df - 1) * PDF *
    ifelse(is.finite(z_lb_ub), 1 + z_lb_ub^2/df, 0)
  CRPS <- crps.t(z, df, 0, 1, lb, ub, "trunc", "trunc")
  
  Hessian <- calcHess_trunc(z, scale, lb, ub,
                            CDF, PDF, PDFp_over_PDF, G, CRPS)
  
  return(Hessian)
}

### Hessian (location, scale) for censored logistic distribution
hessclogis <- function(y, location, scale,
                      lower = -Inf, upper = Inf) {
  
  z <- (pmin(pmax(y, lower), upper) - location) / scale
  lb <- (lower - location) / scale
  ub <- (upper - location) / scale
  
  term1 <- dlogis(z)
  term2 <- dlogis(lb) * plogis(lb)
  term3 <- dlogis(ub) * plogis(ub, lower.tail = FALSE)
  
  d2mu <- term1 - term2 - term3
  dmu.dsigma <- dsigma.dmu <- term1 * z -
    term2 * ifelse(is.finite(lb), lb, 0) -
    term3 * ifelse(is.finite(ub), ub, 0)
  d2sigma <- term1 * z^2 -
    term2 * ifelse(is.finite(lb), lb^2, 0) -
    term3 * ifelse(is.finite(ub), ub^2, 0)
  
  return(2/scale * cbind(d2mu, d2sigma, dmu.dsigma, dsigma.dmu))
}

### Hessian (location, scale) for truncated logistic distribution
hesstlogis <- function(y, location, scale,
                       lower = -Inf, upper = Inf) {
  
  z <- (pmin(pmax(y, lower), upper) - location) / scale
  lb <- (lower - location) / scale
  ub <- (upper - location) / scale
  z_lb_ub <- cbind(z, lb, ub)
  
  CDF <- plogis(z_lb_ub)
  denom <- CDF[, "ub"] - CDF[, "lb"]
  CDF <- CDF / denom
  PDF <- dlogis(z_lb_ub) / denom
  PDFp_over_PDF <- - 1 + 2 * (1 + exp(-z_lb_ub)) * dlogis(z_lb_ub)
  G <- ifelse(is.finite(z_lb_ub),
              plogis(z_lb_ub, log = TRUE) -
                z_lb_ub * plogis(z_lb_ub, lower.tail = FALSE),
              0) / denom
  CRPS <- crps.logis(z, 0, 1, lb, ub, "trunc", "trunc")
  
  Hessian <- calcHess_trunc(z, scale, lb, ub,
                            CDF, PDF, PDFp_over_PDF, G, CRPS)
  
  return(Hessian)
}


## helper function for truncated distributions
calcHess_trunc <- function(z, scale, lb, ub,
                           CDF, PDF, PDFp_over_PDF, G, CRPS) {
  a <- ifelse(is.finite(lb), lb, 0)
  b <- ifelse(is.finite(ub), ub, 0)
  PDFquot_a <- ifelse(is.finite(lb), PDFp_over_PDF[, "lb"], 0)
  PDFquot_b <- ifelse(is.finite(ub), PDFp_over_PDF[, "ub"], 0)
  
  Fxa <- CDF[, "lb"] - CDF[, "z"]
  Fxb <- CDF[, "ub"] - CDF[, "z"]
  Gxa <- G[, "lb"] - G[, "z"]
  Gxb <- G[, "ub"] - G[, "z"]
  
  T_y <- CDF[, "lb"] + CDF[, "ub"] - 2 * CDF[, "z"]
  T_a2 <- (z * Fxa - Gxa + CRPS)
  T_a <- -2 * PDF[, "lb"] * T_a2
  T_b2 <- (z * Fxb - Gxb + CRPS)
  T_b <- 2 * PDF[, "ub"] * T_b2
  
  T_yy <- PDF[, "z"]
  T_ya <- -PDF[, "lb"] * Fxb
  T_yb <- PDF[, "ub"] * Fxa
  #T_ay <- T_ya
  T_aa <- PDF[, "lb"]^2 * (z * (Fxb - Fxa) - CRPS + 4 * T_a2 - a) -
    0.5 * PDFquot_a * T_a
  T_ab <- -PDF[, "lb"] * PDF[, "ub"] * (-CRPS + 2 * T_a2 + 2 * T_b2)
  #T_by <- T_yb
  #T_ba <- T_ab
  T_bb <- PDF[, "ub"]^2 * (z * (Fxa - Fxb) - CRPS + 4 * T_b2 + b) -
    0.5 * PDFquot_b * T_b
  
  d2mu <- 2 / scale * (
    T_yy + T_ya + T_yb +
    T_ya + T_aa + T_ab +
    T_yb + T_ab + T_bb
  )
    
  dmu.dsigma <- dsigma.dmu <- 2 / scale * (
    z * (T_yy + T_ya + T_yb) +
    a * (T_ya + T_aa + T_ab) +
    b * (T_yb + T_ab + T_bb)
  )
    
  d2sigma <- 2 / scale * (
    z^2 * T_yy +
    a^2 * T_aa +
    b^2 * T_bb +
    z * a * 2 * T_ya +
    z * b * 2 * T_yb +
    a * b * 2 * T_ab
  )
  
  Hessian <- cbind(d2mu, d2sigma, dmu.dsigma, d2sigma)
  rownames(Hessian) <- NULL
  
  return(Hessian)
}
