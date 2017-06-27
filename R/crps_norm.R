### crps ###

# standard
#' @export
crps_norm <- function(y, mean = 0, sd = 1) {
  if (identical(mean, 0) && identical(sd , 1)) {
    y * (2 * pnorm(y) - 1) + 2 * dnorm(y) - 1 / sqrt(pi)
  } else if (all(is.finite(sd) & sd > 0)) {
    sd * crps_norm((y - mean) / sd)
  } else {
    input <- data.frame(z = abs(y - mean), sd = sd)
    out <- rep(NaN, dim(input)[1L])
    isNaN <- with(input, is.na(sd) | sd < 0)
    ind1 <- input$sd == 0 & !isNaN
    ind2 <- !isNaN & !ind1
    if (any(ind1)) out[ind1] <- input$z[ind1]
    if (any(ind2)) {
      out[ind2] <- with(input[ind2, ],
                        sd * crps_norm(-z / sd))
    }
    out
  }
}


# censored
#' @export
crps_cnorm <- function(y, location = 0, scale = 1,
                       lower = -Inf, upper = Inf) {
  nan_in_bounds <- anyNA(lower) || anyNA(upper)
  
  if (!nan_in_bounds &&
      identical(location, 0) &&
      identical(scale, 1)) {
    
    out_l1 <- out_u1 <- out_l2 <- 0
    out_u2 <- 1
    z <- y
    if (!identical(lower, -Inf)) {
      p_l <- pnorm(lower)
      out_l1 <- -lower * p_l^2 - 2 * dnorm(lower) * p_l
      out_l1[lower == -Inf] <- 0
      out_l2 <- pnorm(lower, sd = sqrt(0.5))
      z <- pmax(lower, z)
    }
    if (!identical(upper, Inf)) {
      p_u <- pnorm(upper, lower.tail = FALSE)
      out_u1 <- upper * p_u^2 - 2 * dnorm(upper) * p_u
      out_u1[upper == Inf] <- 0
      out_u2 <- pnorm(upper, sd = sqrt(0.5))
      z <- pmin(upper, z)
    }
    b <- out_u2 - out_l2
    out_z <- z * (2 * pnorm(z) - 1) + 2 * dnorm(z)
    
    out <- out_z + out_l1 + out_u1 - b / sqrt(pi)
    out[lower > upper] <- NaN
    out[lower == upper] <- 0
    out + abs(y - z)
  } else if (!nan_in_bounds &&
             all(is.finite(scale) & scale > 0)) {
    if (!identical(lower, -Inf)) lower <- (lower - location) / scale
    if (!identical(upper,  Inf)) upper <- (upper - location) / scale
    scale * crps_cnorm((y - location) / scale,
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
        scale * crps_cnorm(z / scale,
                           lower = lower / scale,
                           upper = upper / scale)
      })
    }
    out
  }
}


# truncated
#' @export
crps_tnorm <- function(y, location = 0, scale = 1,
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
      p_l <- pnorm(lower)
      out_l <- pnorm(lower, sd = sqrt(0.5))
      z <- pmax(lower, z)
    }
    if (!identical(upper,  Inf)) {
      p_u <- pnorm(upper)
      out_u <- pnorm(upper, sd = sqrt(0.5))
      z <- pmin(upper, z)
    }
    a <- p_u - p_l
    b <- out_u - out_l
    b[b == 0] <- NaN
    out_z <- z * (2 * pnorm(z) - p_l - p_u) + 2 * dnorm(z)
    
    out <- (out_z - b / a / sqrt(pi)) / a
    out[lower > upper] <- NaN
    out[lower == upper] <- 0
    out + abs(y - z)
  } else if (!nan_in_bounds &&
             all(is.finite(scale) & scale > 0)) {
    if (!identical(lower, -Inf)) lower <- (lower - location) / scale
    if (!identical(upper,  Inf)) upper <- (upper - location) / scale
    scale * crps_tnorm((y - location) / scale,
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
        scale * crps_tnorm(z / scale,
                           lower = lower / scale,
                           upper = upper / scale)
      })
    }
    out
  }
}


### gradient (location, scale) ###

# standard
#' @export
gradcrps_norm <- function(y , mean = 0, sd = 1) {
  if (identical(mean, 0) &&
      identical(sd, 1)) {
    
    term0 <- crps_norm(y)
    term1 <- 1 - 2 * pnorm(y)
    
    cbind(dmean = term1, dsd = term0 + y * term1)
  } else if (all(is.finite(sd) & sd > 0)) {
    gradcrps_norm((y - mean) / sd)
  } else {
    input <- data.frame(z = y - mean, sd = sd)
    out <- matrix(NaN, dim(input)[1L], 2,
                  dimnames = list(NULL, c("dmean", "dsd")))
    isNaN <- with(input, is.na(sd) | sd <= 0)
    ind2 <- !isNaN
    if (any(ind2)) {
      out[ind2] <- with(input[ind2, ],
                        gradcrps_norm(z / sd))
    }
    out
  }
}


# censored
#' @export
gradcrps_cnorm <- function(y, location = 0, scale = 1,
                           lower = -Inf, upper = Inf) {
  nan_in_bounds <- anyNA(lower) || anyNA(upper)
  
  if (!nan_in_bounds &&
      identical(location, 0) &&
      identical(scale, 1)) {
    
    z <- y
    if (!identical(lower, -Inf)) z <- pmax(lower, z)
    if (!identical(upper,  Inf)) z <- pmin(upper, z)
    
    term0 <- crps_cnorm(z, lower = lower, upper = upper)
    term1 <- 1 - 2 * pnorm(z)
    term2 <- pnorm(lower)^2
    term3 <- pnorm(upper, lower.tail = FALSE)^2
    
    dloc <- term1 + term2 - term3
    dscale <- term0 + term1 * z +
      term2 * ifelse(lower == -Inf, lower, 0) -
      term3 * ifelse(upper ==  Inf, upper, 0)
    
    out <- cbind(dloc, dscale)
    out[lower > upper, ] <- NaN
    out[lower == upper, ] <- 0
    out
  } else if (!nan_in_bounds &&
             all(is.finite(scale) & scale > 0)) {
    if (!identical(lower, -Inf)) lower <- (lower - location) / scale
    if (!identical(upper,  Inf)) upper <- (upper - location) / scale
    gradcrps_cnorm((y - location) / scale,
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
        gradcrps_cnorm(z / scale,
                       lower = lower / scale,
                       upper = upper / scale)
      })
    }
    out
  }
}


# truncated
#' @export
gradcrps_tnorm <- function(y, location = 0, scale = 1,
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
    
    p_l <- d_l <- d_u <- 0
    p_u <- 1
    z <- y
    if (!identical(lower, -Inf)) {
      p_l <- pnorm(lower)
      d_l <- dnorm(lower)
      z <- pmax(lower, z)
    }
    if (!identical(upper, Inf)) {
      p_u <- pnorm(upper)
      d_u <- dnorm(upper)
      z <- pmin(upper, z)
    }
    a <- p_u - p_l
    p_z <- pnorm(z)
    d_z <- dnorm(z)
    
    term0 <- crps_tnorm(z, lower = lower, upper = upper)
    term1 <- (p_u + p_l - 2 * p_z) / a
    term2 <- 2 * d_l / a * 
      ((z * p_l - z * p_z + d_l - d_z) / a + term0)
    term3 <- 2 * d_u / a *
      ((z * p_u - z * p_z + d_u - d_z) / a + term0)
    
    dloc <- (term1 - term2 + term3) * sign
    dscale <- term0 + term1 * z - 
      term2 * ifelse(lower == -Inf, lower, 0) +
      term3 * ifelse(upper ==  Inf, upper, 0)
    
    out <- cbind(dloc, dscale)
    out[lower > upper, ] <- NaN
    out[lower == upper, ] <- 0
    out
  } else if (!nan_in_bounds &&
             all(is.finite(scale) & scale > 0)) {
    if (!identical(lower, -Inf)) lower <- (lower - location) / scale
    if (!identical(upper,  Inf)) upper <- (upper - location) / scale
    gradcrps_tnorm((y - location) / scale,
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
        gradcrps_tnorm(z / scale,
                       lower = lower / scale,
                       upper = upper / scale)
      })
    }
    out
  }
}


### Hessian (location, scale) ###

# standard
#' @export
hesscrps_norm <- function(y , mean = 0, sd = 1) {
  if (identical(mean, 0) &&
      identical(sd, 1)) {
    
    term1 <- dnorm(y)
    
    d2mean <- term1
    dmean.dsd <- dsd.dmean <- term1 * y
    d2sd <- term1 * y^2
    
    2 * cbind(d2mean, d2sd, dmean.dsd, dsd.dmean)
  } else if (all(is.finite(sd) & sd > 0)) {
    hesscrps_norm((y - mean) / sd) / sd
  } else {
    input <- data.frame(z = y - mean, sd = sd)
    out <- matrix(NaN, dim(input)[1L], 4,
                  dimnames = list(NULL, c("d2mean", "d2sd",
                                          "dmean.dsd", "dsd.dmean")))
    isNaN <- with(input, is.na(sd) | sd <= 0)
    ind2 <- !isNaN
    if (any(ind2)) {
      out[ind2] <- with(input[ind2, ],
                        hesscrps_norm(z / sd) / sd)
    }
    out
  }
}


# censored
#' @export
hesscrps_cnorm <- function(y, location = 0, scale = 1,
                           lower = -Inf, upper = Inf) {
  nan_in_bounds <- anyNA(lower) || anyNA(upper)
  
  if (!nan_in_bounds &&
      identical(location, 0) &&
      identical(scale, 1)) {
    
    z <- y
    if (!identical(lower, -Inf)) z <- pmax(lower, z)
    if (!identical(upper,  Inf)) z <- pmin(upper, z)
    
    term1 <- dnorm(z)
    term2 <- dnorm(lower) * pnorm(lower)
    term3 <- dnorm(upper) * pnorm(upper, lower.tail = FALSE)
    
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
    hesscrps_cnorm((y - location) / scale,
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
        hesscrps_cnorm(z / scale,
                       lower = lower / scale,
                       upper = upper / scale) / scale
      })
    }
    out
  }
}
  

# truncated
#' @export
hesscrps_tnorm <- function(y, location = 0, scale = 1,
                      lower = -Inf, upper = Inf) {
  nan_in_bounds <- anyNA(lower) || anyNA(upper)
  
  if (!nan_in_bounds &&
      identical(location, 0) &&
      identical(scale, 1)) {
    
    z <- y
    if (!identical(lower, -Inf)) z <- pmax(lower, z)
    if (!identical(upper,  Inf)) z <- pmin(upper, z)
    
    z_lb_ub <- cbind(z, lb = lower, ub = upper)
    
    CDF <- pnorm(z_lb_ub)
    denom <- CDF[, "ub"] - CDF[, "lb"]
    CDF <- CDF / denom
    PDF <- dnorm(z_lb_ub) / denom
    PDFp_over_PDF <- -z_lb_ub
    G <- -PDF
    CRPS <- crps_tnorm(z, lower = lower, upper = upper)
    
    out <- calcHess_trunc(z, 1, lower, upper,
                          CDF, PDF, PDFp_over_PDF, G, CRPS)
    out[lower > upper, ] <- NaN
    out[lower == upper, ] <- 0
    out
  } else if (!nan_in_bounds &&
             all(is.finite(scale) & scale > 0)) {
    if (!identical(lower, -Inf)) lower <- (lower - location) / scale
    if (!identical(upper,  Inf)) upper <- (upper - location) / scale
    hesscrps_tnorm((y - location) / scale,
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
        hesscrps_cnorm(z / scale,
                       lower = lower / scale,
                       upper = upper / scale) / scale
      })
    }
    out
  }
}