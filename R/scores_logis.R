#' Calculating scores for the logistic distribution
#' 
#' These functions calculate scores (CRPS, logarithmic score) and its gradient and Hessian with respect
#' to the parameters of a location-scale transformed logistic
#' distribution. Furthermore, the censoring transformation and
#' the truncation transformation may be introduced on top of the
#' location-scale transformed logistic distribution.
#' 
#' @usage
#' ## score functions
#' crps_logis(y, location = 0, scale = 1)
#' crps_clogis(y, location = 0, scale = 1, lower = -Inf, upper = Inf)
#' crps_tlogis(y, location = 0, scale = 1, lower = -Inf, upper = Inf)
#' crps_gtclogis(y, location = 0, scale = 1, lower = -Inf, upper = Inf, lmass = 0, umass = 0)
#' logs_logis(y, location = 0, scale = 1)
#' logs_tlogis(y, location = 0, scale = 1, lower = -Inf, upper = Inf)
#' dss_logis(y, location = 0, scale = 1)
#'
#' ## gradient (location, scale) functions
#' gradcrps_logis(y, location = 0, scale = 1)
#' gradcrps_clogis(y, location = 0, scale = 1, lower = -Inf, upper = Inf)
#' gradcrps_tlogis(y, location = 0, scale = 1, lower = -Inf, upper = Inf)
#'
#' ## Hessian (location, scale) functions
#' hesscrps_logis(y, location = 0, scale = 1)
#' hesscrps_clogis(y, location = 0, scale = 1, lower = -Inf, upper = Inf)
#' hesscrps_tlogis(y, location = 0, scale = 1, lower = -Inf, upper = Inf)
#' 
#' @param y vector of observations.
#' @param location vector of location parameters.
#' @param scale vector of scale paramters.
#' @param lower,upper lower and upper truncation/censoring bounds.
#' @param lmass,umass vectors of point masses in \code{lower} and \code{upper}
#'  respectively. 
#' @return For the score functions: a vector of score values.
#' 
#' For the gradient and Hessian functions: a matrix with column names
#' corresponding to the respective partial derivatives.
#' @name scores_logis
#' @importFrom stats plogis
NULL

# logistic distribution
# z = plogis(z, log.p = TRUE) - plogis(-z, log.p = TRUE)
# Gz = z * plogis(z) + plogis(-z, log.p = TRUE)
#
# Taylor series expansion for z -> -Inf
# plogis(z) + plogis(-z, log.p = TRUE) = -plogis(z)^2/2 - plogis(z)^3/3

### crps ###

# standard
#' @rdname scores_logis
#' @usage NULL
#' @export
crps_logis <- function(y, location = 0, scale = 1) {
  if (identical(location, 0) & identical(scale, 1)) {
    y - 2 * plogis(y, log.p = TRUE) - 1
  } else if (all(is.finite(scale) & scale > 0)) {
    scale * crps_logis((y - location) / scale)
  } else {
    input <- data.frame(z = abs(y - location), scale = scale)
    out <- rep(NaN, dim(input)[1L])
    isNaN <- with(input, is.na(scale) | scale < 0)
    ind1 <- !isNaN & input$scale == 0
    ind2 <- !isNaN & !ind1
    if (any(ind1)) out[ind1] <- input$z[ind1]
    if (any(ind2)) {
      out[ind2] <- with(input[ind2, ],
                        scale * crps_logis(-z / scale))
    }
    out
  }
}


# censored
#' @rdname scores_logis
#' @usage NULL
#' @export
crps_clogis <- function(y, location = 0, scale = 1,
                       lower = -Inf, upper = Inf) {
  
  nan_in_bounds <- anyNA(lower) || anyNA(upper)
  
  if (!nan_in_bounds &&
      identical(location, 0) &&
      identical(scale, 1)) {
    
    out_l <- out_u <- 0
    z <- y
    if (!identical(lower, -Inf)) {
      p_l <- plogis(lower)
      lp_ml <- plogis(-lower, log.p = TRUE)
      out_l <- p_l + lp_ml
      # Taylor series expansion of 'x + log(1-x)' at 0 to avoid underflow
      # out_l <- ifelse(p_l > 1e-8, p_l + lp_ml, -p_l^2/2 - p_l^3/3)
      z <- pmax(lower, z)
    }
    if (!identical(upper, Inf)) {
      p_mu <- plogis(-upper)
      lp_u <- plogis(upper, log.p = TRUE)
      out_u <- p_mu + lp_u
      # out_u <- ifelse(p_mu > 1e-8, p_mu + lp_u, -p_mu^2/2 - p_mu^3/3)
      z <- pmin(upper, z)
    }
    out_z <- z - 2 * plogis(z, log.p = TRUE) - 1
    
    out <- out_z + out_l + out_u
    out[lower > upper] <- NaN
    out[lower == upper] <- 0
    out + abs(y - z)
  } else if (!nan_in_bounds &&
             all(is.finite(scale) & scale > 0)) {
    if (!identical(lower, -Inf)) lower <- (lower - location) / scale
    if (!identical(upper,  Inf)) upper <- (upper - location) / scale
    scale * crps_clogis((y - location) / scale,
                        lower = lower, upper = upper)
  } else {
    input <- data.frame(z = y - location, scale = scale,
                        lower = lower - location,
                        upper = upper - location)
    out <- rep(NaN, dim(input)[1L])
    isNaN <- with(input, {
      is.na(scale) | scale < 0 |
        is.na(lower) | is.na(upper)
    })
    ind1 <- !isNaN & with(input, scale == 0 & lower <= upper)
    ind2 <- !isNaN & input$scale > 0
    if (any(ind1)) {
      out[ind1] <- with(input[ind1, ], {
        abs(z - pmax(lower, 0) - pmin(upper, 0))
      })
    }
    if (any(ind2)) {
      out[ind2] <- with(input[ind2, ], {
        scale * crps_clogis(z / scale,
                           lower = lower / scale,
                           upper = upper / scale)
      })
    }
    out
  }
}


# truncated
#' @rdname scores_logis
#' @usage NULL
#' @export
crps_tlogis <- function(y, location = 0, scale = 1,
                        lower = -Inf, upper = Inf) {
  nan_in_bounds <- anyNA(lower) || anyNA(upper)
  
  if (!nan_in_bounds &&
      identical(location, 0) &&
      identical(scale, 1)) {
    
    ind_swap <- lower > 3
    if (any(ind_swap)) {
      sign <- 1 - 2 * ind_swap
      y <- y * sign
      lower <- lower * sign
      upper <- upper * sign
      l_tmp <- lower[ind_swap]
      u_tmp <- upper[ind_swap]
      lower[ind_swap] <- u_tmp
      upper[ind_swap] <- l_tmp
    }
    
    out_l <- p_l <- 0
    out_u <- p_u <- 1
    z <- y
    if (!identical(lower, -Inf)) {
      p_l <- plogis(lower)
      lp_ml <- plogis(-lower, log.p = TRUE)
      # Taylor series expansion of 'x + log(1-x)' at 0 to avoid underflow
      out_l <- ifelse(p_l > 1e-8, p_l + lp_ml, -p_l^2/2 - p_l^3/3) -
        p_l * (lower * p_l + 2 * lp_ml)
      out_l[lower == -Inf] <- 0
      z <- pmax(lower, z)
    }
    if (!identical(upper, Inf)) {
      p_u <- plogis(upper)
      lp_mu <- plogis(-upper, log.p = TRUE)
      out_u <- ifelse(p_u > 1e-8, p_u + lp_mu, -p_u^2/2 - p_u^3/3) -
        p_u * (upper * p_u + 2 * lp_mu)
      out_u[upper ==  Inf] <- 1
      z <- pmin(upper, z)
    }
    a <- p_u - p_l
    b <- out_u - out_l
    b[b == 0] <- NaN
    out_z <- z * (-p_l - p_u) - 2 * plogis(-z, log.p = TRUE)
    
    out <- (out_z - b / a) / a
    out[lower > upper] <- NaN
    out[lower == upper] <- 0
    out + abs(y - z)
  } else if (!nan_in_bounds &&
             all(is.finite(scale) & scale > 0)) {
    if (!identical(lower, -Inf)) lower <- (lower - location) / scale
    if (!identical(upper,  Inf)) upper <- (upper - location) / scale
    scale * crps_tlogis((y - location) / scale,
                        lower = lower, upper = upper)
  } else {
    input <- data.frame(z = y - location, scale = scale,
                        lower = lower - location,
                        upper = upper - location)
    out <- rep(NaN, dim(input)[1L])
    isNaN <- with(input, {
      is.na(scale) | scale < 0 |
        is.na(lower) | is.na(upper)
    })
    ind1 <- !isNaN & with(input, scale == 0 & lower < 0 & upper > 0)
    ind2 <- !isNaN & input$scale > 0
    if (any(ind1)) out[ind1] <- abs(input$z[ind1])
    if (any(ind2)) {
      out[ind2] <- with(input[ind2, ], {
        scale * crps_tlogis(z / scale,
                            lower = lower / scale,
                            upper = upper / scale)
      })
    }
    out
  }
}


# generalized truncated/censored
#' @rdname scores_logis
#' @usage NULL
#' @export
crps_gtclogis <- function(y, location = 0, scale = 1,
                          lower = -Inf, upper = Inf,
                          lmass = 0, umass = 0) {
  if (!identical(location, 0)) {
    y <- y - location
    if (!identical(lower, -Inf)) lower <- lower - location
    if (!identical(upper,  Inf)) upper <- upper - location
  }
  
  if (identical(scale, 1)) {
    ind_swap <- lower > 3
    if (any(ind_swap, na.rm = TRUE)) {
      sign <- 1 - 2 * ind_swap
      y <- y * sign
      lower <- lower * sign
      upper <- upper * sign
      l_tmp <- lower[ind_swap]
      u_tmp <- upper[ind_swap]
      lower[ind_swap] <- u_tmp
      upper[ind_swap] <- l_tmp
      if (length(lmass) < length(lower))
        lmass <- rep_len(lmass, length(lower))
      if (length(umass) < length(upper))
        umass <- rep_len(umass, length(upper))
      l_tmp <- lmass[ind_swap]
      u_tmp <- umass[ind_swap]
      lmass[ind_swap] <- u_tmp
      umass[ind_swap] <- l_tmp
    }
    
    out_l1 <- out_l2 <- out_l3 <- out_u1 <- out_u2 <- p_l <- 0
    out_u3 <- p_u <- 1
    z <- y
    if (!identical(lower, -Inf) || !identical(lmass, 0)) {
      lmass[lmass < 0 & lmass > 1] <- NaN
      p_l <- plogis(lower)
      lp_ml <- plogis(-lower, log.p = TRUE)
      ind <- lower == -Inf
      out_l1 <- lower * lmass^2
      out_l1[lmass == 0] <- 0
      out_l2 <- 2 * (lower * p_l + lp_ml) * lmass
      out_l2[ind] <- 0
      # Taylor series expansion of 'x + log(1-x)' at 0 to avoid underflow
      out_l3 <- ifelse(p_l > 1e-8, p_l + lp_ml, -p_l^2/2 - p_l^3/3) -
        p_l * (lower * p_l + 2 * lp_ml)
      out_l3[ind] <- 0
      z <- pmax(lower, z)
    }
    if (!identical(upper, Inf) || !identical(umass, 0)) {
      umass[umass < 0 & umass > 1] <- NaN
      p_u <- plogis(upper)
      lp_mu <- plogis(-upper, log.p = TRUE)
      ind <- upper == Inf
      out_u1 <- upper * umass^2
      out_u1[umass == 0] <- 0
      out_u2 <- 2 * (upper * p_u + lp_mu) * umass
      out_u2[ind] <- 0
      out_u3 <- ifelse(p_u > 1e-8, p_u + lp_mu, -p_u^2/2 - p_u^3/3) -
        p_u * (upper * p_u + 2 * lp_mu)
      out_u3[ind] <- 1
      z <- pmin(upper, z)
    }
    
    a1 <- p_u - p_l
    a2 <- 1 - (umass + lmass)
    a2[a2 < 0 | a2 > 1] <- NaN
    b <- out_u3 - out_l3
    b[b == 0] <- NaN
    
    out <- out_u1 - out_l1 -
      (z * ((1 - 2 * lmass) * p_u + (1 - 2 * umass) * p_l) +
         (2 * plogis(-z, log.p = TRUE) - out_u2 - out_l2 + a2 * b / a1) * a2
      ) / a1
    
    out[lower > upper] <- NaN
    out[lower == upper] <- 0
    out + abs(y - z)
  } else {
    scale[scale < 0] <- NaN
    if (!identical(lower, -Inf)) lower <- lower / scale
    if (!identical(upper,  Inf)) upper <- upper / scale
    if (all(scale > 0, na.rm = TRUE)) {
      scale * crps_gtclogis(y / scale,
                            lower = lower, upper = upper,
                            lmass = lmass, umass = umass)
    } else {
      out <- scale * crps_gtclogis(y / scale,
                                   lower = lower, upper = upper,
                                   lmass = lmass, umass = umass)
      ind <- scale == 0 & lower < 0 & upper > 0
      out[ind] <- rep_len(
        (pmin(y, 0) - lower) * lmass^2 - pmin(y, 0) * (1 - lmass)^2 +
          (upper - pmax(y, 0)) * umass^2 + pmax(y, 0) * (1 - umass)^2,
        length(out)
      )[ind]
      out
    }
  }
}

#' @rdname scores_logis
#' @usage NULL
#' @export
logs_logis <- function(y, location = 0, scale = 1) {
  -flogis(y, location, scale, log = TRUE)
}

#' @rdname scores_logis
#' @usage NULL
#' @export
logs_tlogis <- function(y, location = 0, scale = 1,
                        lower = -Inf, upper = Inf) {
  -flogis(y, location, scale, lower, upper, log = TRUE)
}

#' @rdname scores_logis
#' @usage NULL
#' @export
dss_logis <- function(y, location = 0, scale = 1) {
  if (!identical(location, 0)) y <- y - location
  scale[scale <= 0] <- NaN
  s <- scale * pi / sqrt(3)
  (y / s)^2 + 2 * log(s)
}


### gradients (location, scale)

# standard
#' @rdname scores_logis
#' @usage NULL
#' @export
gradcrps_logis <- function(y , location = 0, scale = 1) {
  if (identical(location, 0) &&
      identical(scale, 1)) {
    
    term0 <- crps_logis(y)
    term1 <- 1 - 2 * plogis(y)
    
    cbind(dloc = term1, dscale = term0 + y * term1)
  } else if (all(is.finite(scale) & scale > 0)) {
    gradcrps_logis((y - location) / scale)
  } else {
    input <- data.frame(z = y - location, scale = scale)
    out <- matrix(NaN, dim(input)[1L], 2,
                  dimnames = list(NULL, c("dloc", "dscale")))
    isNaN <- with(input, is.na(scale) | scale <= 0)
    ind2 <- !isNaN
    if (any(ind2)) {
      out[ind2, ] <- with(input[ind2, ], {
        gradcrps_logis(z / scale)
      })
    }
    out
  }
}

# censored
#' @rdname scores_logis
#' @usage NULL
#' @export
gradcrps_clogis <- function(y, location = 0, scale = 1,
                            lower = -Inf, upper = Inf) {
  nan_in_bounds <- anyNA(lower) || anyNA(upper)
  
  if (!nan_in_bounds &&
      identical(location, 0) &&
      identical(scale, 1)) {
    
    z <- y
    if (!identical(lower, -Inf)) z <- pmax(lower, z)
    if (!identical(upper,  Inf)) z <- pmin(upper, z)
    
    term0 <- crps_clogis(z, lower = lower, upper = upper)
    term1 <- 1 - 2 * plogis(z)
    term2 <- plogis(lower)^2
    term3 <- plogis(upper, lower.tail = FALSE)^2
    
    dloc <- term1 + term2 - term3
    dscale <- term0 + term1 * z +
      term2 * ifelse(lower == -Inf, 0, lower) -
      term3 * ifelse(upper ==  Inf, 0, upper)
    
    out <- cbind(dloc, dscale)
    out[lower > upper, ] <- NaN
    out[lower == upper, ] <- 0
    out
  } else if (!nan_in_bounds &&
             all(is.finite(scale) & scale > 0)) {
    if (!identical(lower, -Inf)) lower <- (lower - location) / scale
    if (!identical(upper,  Inf)) upper <- (upper - location) / scale
    gradcrps_clogis((y - location) / scale,
                    lower = lower, upper = upper)
  } else {
    input <- data.frame(z = y - location, scale = scale,
                        lower = lower - location,
                        upper = upper - location)
    out <- matrix(NaN, dim(input)[1L], 2,
                  dimnames = list(NULL, c("dloc", "dscale")))
    isNaN <- with(input, {
      is.na(scale) | scale <= 0 |
        is.na(lower) | is.na(upper)
    })
    ind2 <- !isNaN
    if (any(ind2)) {
      out[ind2, ] <- with(input[ind2, ], {
        gradcrps_clogis(z / scale,
                       lower = lower / scale,
                       upper = upper / scale)
      })
    }
    out
  }
}


# truncated
#' @rdname scores_logis
#' @usage NULL
#' @export
gradcrps_tlogis <- function(y, location = 0, scale = 1,
                            lower = -Inf, upper = Inf) {
  nan_in_bounds <- anyNA(lower) || anyNA(upper)
  
  if (!nan_in_bounds &&
      identical(location, 0) &&
      identical(scale, 1)) {
    
    sign <- 1
    ind_swap <- lower > 3
    if (any(ind_swap)) {
      sign <- 1 - 2 * ind_swap
      y <- y * sign
      lower <- lower * sign
      upper <- upper * sign
      l_tmp <- lower[ind_swap]
      u_tmp <- upper[ind_swap]
      lower[ind_swap] <- u_tmp
      upper[ind_swap] <- l_tmp
    }
    
    p_l <- d_l <- d_u <- G_l <- G_u <- 0
    p_u <- 1
    z <- y
    if (!identical(lower, -Inf)) {
      p_l <- plogis(lower)
      d_l <- dlogis(lower)
      G_l <- lower * p_l + plogis(-lower, log.p = TRUE)
      G_l[is.infinite(lower)] <- 0
      z <- pmax(lower, z)
    }
    if (!identical(upper, Inf)) {
      p_u <- plogis(upper)
      d_u <- dlogis(upper)
      G_u <- upper * p_u + plogis(-upper, log.p = TRUE)
      G_u[is.infinite(upper)] <- 0
      z <- pmin(upper, z)
    }
    a <- p_u - p_l
    p_z <- plogis(z)
    d_z <- dlogis(z)
    Gz_minus_zFz <- plogis(-z, log.p = TRUE)
    
    term0 <- crps_tlogis(z, lower = lower, upper = upper)
    term1 <- (p_l + p_u - 2 * p_z) / a
    term2 <- 2 * d_l / a * 
      ((z * p_l - G_l + Gz_minus_zFz) / a + term0)
    term3 <- 2 * d_u / a *
      ((z * p_u - G_u + Gz_minus_zFz) / a + term0)
    
    dloc <- (term1 - term2 + term3) * sign
    dscale <- term0 + term1 * z - 
      term2 * ifelse(lower == -Inf, 0, lower) +
      term3 * ifelse(upper ==  Inf, 0, upper)
    
    out <- cbind(dloc, dscale)
    out[lower > upper, ] <- NaN
    out[lower == upper, ] <- 0
    out
  } else if (!nan_in_bounds &&
             all(is.finite(scale) & scale > 0)) {
    if (!identical(lower, -Inf)) lower <- (lower - location) / scale
    if (!identical(upper,  Inf)) upper <- (upper - location) / scale
    gradcrps_tlogis((y - location) / scale,
                    lower = lower, upper = upper)
  } else {
    input <- data.frame(z = y - location, scale = scale,
                        lower = lower - location,
                        upper = upper - location)
    out <- matrix(NaN, dim(input)[1L], 2,
                  dimnames = list(NULL, c("dloc", "dscale")))
    isNaN <- with(input, {
      is.na(scale) | scale <= 0 |
        is.na(lower) | is.na(upper)
    })
    ind2 <- !isNaN
    if (any(ind2)) {
      out[ind2, ] <- with(input[ind2, ], {
        gradcrps_tlogis(z / scale,
                        lower = lower / scale,
                        upper = upper / scale)
      })
    }
    out
  }
}


### Hessian (location, scale) ###

# standard
#' @rdname scores_logis
#' @usage NULL
#' @export
hesscrps_logis <- function(y , location = 0, scale = 1) {
  if (identical(location, 0) &&
      identical(scale, 1)) {
    
    term1 <- dlogis(y)
    
    d2mu <- term1
    dmu.dsigma <- dsigma.dmu <- term1 * y
    d2sigma <- term1 * y^2
    
    2 * cbind(d2mu, d2sigma, dmu.dsigma, dsigma.dmu)
  } else if (all(is.finite(scale) & scale > 0)) {
    hesscrps_logis((y - location) / scale) / scale
  } else {
    input <- data.frame(z = y - location, scale = scale)
    out <- matrix(NaN, dim(input)[1L], 4,
                  dimnames = list(NULL, c("d2loc", "d2scale",
                                          "dloc.dscale", "dscale.dloc")))
    isNaN <- with(input, is.na(scale) | scale <= 0)
    ind2 <- !isNaN
    if (any(ind2)) {
      out[ind2, ] <- with(input[ind2, ], {
        hesscrps_logis(z / scale) / scale
      })
    }
    out
  }
}

# censored
#' @rdname scores_logis
#' @usage NULL
#' @export
hesscrps_clogis <- function(y, location = 0, scale = 1,
                            lower = -Inf, upper = Inf) {
  nan_in_bounds <- anyNA(lower) || anyNA(upper)
  
  if (!nan_in_bounds &&
      identical(location, 0) &&
      identical(scale, 1)) {
    
    z <- y
    if (!identical(lower, -Inf)) z <- pmax(lower, z)
    if (!identical(upper,  Inf)) z <- pmin(upper, z)
    
    term1 <- dlogis(z)
    term2 <- dlogis(lower) * plogis(lower)
    term3 <- dlogis(upper) * plogis(upper, lower.tail = FALSE)
    
    d2loc <- term1 - term2 - term3
    dloc.dscale <- dscale.dloc <- term1 * z -
      term2 * ifelse(is.finite(lower), lower, 0) -
      term3 * ifelse(is.finite(upper), upper, 0)
    d2scale <- term1 * z^2 -
      term2 * ifelse(is.finite(lower), lower^2, 0) -
      term3 * ifelse(is.finite(upper), upper^2, 0)
    
    out <- 2 * cbind(d2loc, d2scale, dloc.dscale, dscale.dloc)
    out[lower > upper, ] <- NaN
    out[lower == upper, ] <- 0
    out
  } else if (!nan_in_bounds &&
             all(is.finite(scale) & scale > 0)) {
    if (!identical(lower, -Inf)) lower <- (lower - location) / scale
    if (!identical(upper,  Inf)) upper <- (upper - location) / scale
    hesscrps_clogis((y - location) / scale,
                    lower = lower, upper = upper) / scale
  } else {
    input <- data.frame(z = y - location, scale = scale,
                        lower = lower - location,
                        upper = upper - location)
    out <- rep(NaN, dim(input)[1L], 4,
               dimnames = list(NULL, c("d2loc", "d2scale",
                                       "dloc.dscale", "dscale.dloc")))
    isNaN <- with(input, {
      is.na(scale) | scale <= 0 |
        is.na(lower) | is.na(upper)
    })
    ind2 <- !isNaN
    if (any(ind2)) {
      out[ind2] <- with(input[ind2, ], {
        hesscrps_clogis(z / scale,
                        lower = lower / scale,
                        upper = upper / scale) / scale
      })
    }
    out
  }
}


# truncated
#' @rdname scores_logis
#' @usage NULL
#' @export
hesscrps_tlogis <- function(y, location = 0, scale = 1,
                            lower = -Inf, upper = Inf) {
  nan_in_bounds <- anyNA(lower) || anyNA(upper)
  
  if (!nan_in_bounds &&
      identical(location, 0) &&
      identical(scale, 1)) {
    
    z <- y
    if (!identical(lower, -Inf)) z <- pmax(lower, z)
    if (!identical(upper,  Inf)) z <- pmin(upper, z)
    
    z_lb_ub <- cbind(z, lb = lower, ub = upper)
    
    CDF <- plogis(z_lb_ub)
    denom <- CDF[, "ub"] - CDF[, "lb"]
    CDF <- CDF / denom
    PDF <- dlogis(z_lb_ub) / denom
    PDFp_over_PDF <- -1 + 2 * plogis(-z_lb_ub)
    G <- (plogis(-z_lb_ub, log.p = TRUE) +
            z_lb_ub * plogis(z_lb_ub)) / denom
    G[is.infinite(z_lb_ub)] <- 0
    CRPS <- crps_tlogis(z, lower = lower, upper = upper)
    
    out <- calcHess_trunc(z, 1, lower, upper,
                          CDF, PDF, PDFp_over_PDF, G, CRPS)
    out[lower > upper, ] <- NaN
    out[lower == upper, ] <- 0
    out
  } else if (!nan_in_bounds &&
             all(is.finite(scale) & scale > 0)) {
    if (!identical(lower, -Inf)) lower <- (lower - location) / scale
    if (!identical(upper,  Inf)) upper <- (upper - location) / scale
    hesscrps_tlogis((y - location) / scale,
                    lower = lower, upper = upper) / scale
  } else {
    input <- data.frame(z = y - location, scale = scale,
                        lower = lower - location,
                        upper = upper - location)
    out <- rep(NaN, dim(input)[1L], 4,
               dimnames = list(NULL, c("d2loc", "d2scale",
                                       "dloc.dscale", "dscale.dloc")))
    isNaN <- with(input, {
      is.na(scale) | scale <= 0 |
        is.na(lower) | is.na(upper)
    })
    ind2 <- !isNaN
    if (any(ind2)) {
      out[ind2] <- with(input[ind2, ], {
        hesscrps_tlogis(z / scale,
                       lower = lower / scale,
                       upper = upper / scale) / scale
      })
    }
    out
  }
}


################################## Checks ######################################
check_crps_logis <- function(input) {
  required <- c("y", "location", "scale")
  checkNames1(required, names(input))
  checkNumeric(input)
  checkVector(input)
  
  if (any(input$scale <= 0))
    stop("Parameter 'scale' contains non-positive values.")
}

check_crps_clogis <- function(input) {
  required <- c("y", "location", "scale", "lower", "upper")
  checkNames1(required, names(input))
  checkNumeric(input, infinite_exception = c("lower", "upper"))
  checkVector(input)
  
  if (any(input$scale <= 0))
    stop("Parameter 'scale' contains non-positive values.")
  if (any(input$lower > input$upper))
    stop("Parameter 'lower' contains values greater than corresponding values in 'upper'.")
}

check_crps_tlogis <- check_crps_clogis

check_crps_gtclogis <- function(input) {
  required <- c("y", "location", "scale", "lower", "upper", "lmass", "umass")
  checkNames1(required, names(input))
  checkNumeric(input, infinite_exception = c("lower", "upper"))
  checkVector(input)
  
  if (any(input$scale <= 0))
    stop("Parameter 'scale' contains non-positive values.")
  if (any(input$lower > input$upper))
    stop("Parameter 'lower' contains values greater than corresponding values in 'upper'.")
  if (any(input$lmass < 0 | input$lmass > 1))
    stop("Parameter 'lmass' contains values not in [0, 1].")
  if (any(input$umass < 0 | input$umass > 1))
    stop("Parameter 'umass' contains values not in [0, 1].")
  if (any(input$umass + input$lmass > 1))
    stop("Values in 'lmass' and 'umass' add up to more than 1.")
}

check_logs_logis <- check_crps_logis

check_logs_tlogis <- check_crps_tlogis
