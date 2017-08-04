#' Calculating scores for the normal distribution
#' 
#' These functions calculate scores (CRPS, logarithmic score) and their gradient and Hessian with respect
#' to the parameters of a location-scale transformed normal
#' distribution. Furthermore, the censoring transformation and
#' the truncation transformation may be introduced on top of the
#' location-scale transformed normal distribution.
#' 
#' @usage
#' ## score functions
#' crps_norm(y, mean = 0, sd = 1, location = mean, scale = sd)
#' crps_cnorm(y, location = 0, scale = 1, lower = -Inf, upper = Inf)
#' crps_tnorm(y, location = 0, scale = 1, lower = -Inf, upper = Inf)
#' crps_gtcnorm(y, location = 0, scale = 1, lower = -Inf, upper = Inf, lmass = 0, umass = 0)
#' logs_norm(y, mean = 0, sd = 1, location = mean, scale = sd)
#' logs_tnorm(y, location = 0, scale = 1, lower = -Inf, upper = Inf)
#'
#' ## gradient (location, scale) functions
#' gradcrps_norm(y, location = 0, scale = 1)
#' gradcrps_cnorm(y, location = 0, scale = 1, lower = -Inf, upper = Inf)
#' gradcrps_tnorm(y, location = 0, scale = 1, lower = -Inf, upper = Inf)
#'
#' ## Hessian (location, scale) functions
#' hesscrps_norm(y, location = 0, scale = 1)
#' hesscrps_cnorm(y, location = 0, scale = 1, lower = -Inf, upper = Inf)
#' hesscrps_tnorm(y, location = 0, scale = 1, lower = -Inf, upper = Inf)
#' 
#' @param y vector of observations.
#' @param mean an alternative way to specify \code{location}.
#' @param sd an alternative way to specify \code{scale}.
#' @param location vector of location parameters.
#' @param scale vector of scale parameters.
#' @param lower,upper lower and upper truncation/censoring bounds.
#' @param lmass,umass vectors of point masses in \code{lower} and \code{upper}
#'  respectively. 
#' @return For the score functions: a vector of score values.
#' 
#' For the gradient and Hessian functions: a matrix with column names
#' corresponding to the respective partial derivatives.
#' @name scores_norm
#' @importFrom stats pnorm dnorm
NULL


### crps ###

# standard
#' @rdname scores_norm
#' @usage NULL
#' @export
crps_norm <- function(y, mean = 0, sd = 1, location = mean, scale = sd) {
  if (!missing(mean) && !missing(location))
    stop("specify 'mean' or 'location' but not both")
  if (!missing(sd) && !missing(scale))
    stop("specify 'sd' or 'scale' but not both")
  if (!identical(location, 0)) y <- y - location
  if (identical(scale, 1)) {
    y * (2 * pnorm(y) - 1) + (sqrt(2) * exp(-0.5 * y^2) - 1) / sqrt(pi)
  } else {
    z <- y / scale
    z[y == 0 & scale == 0] <- 0
    y * (2 * pnorm(y, sd = scale) - 1) +
      scale * (sqrt(2) * exp(-0.5 * z^2) - 1) / sqrt(pi)
  }
}


# censored
#' @rdname scores_norm
#' @usage NULL
#' @export
crps_cnorm <- function(y, location = 0, scale = 1,
                       lower = -Inf, upper = Inf) {
  if (!identical(location, 0)) {
    y <- y - location
    if (!identical(lower, -Inf)) lower <- lower - location
    if (!identical(upper,  Inf)) upper <- upper - location
  }
  
  if (identical(scale, 1)) {
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
  } else {
    scale[scale < 0] <- NaN
    if (!identical(lower, -Inf)) lower <- lower / scale
    if (!identical(upper,  Inf)) upper <- upper / scale
    if (all(scale > 0, na.rm = TRUE)) {
      scale * crps_cnorm(y / scale, lower = lower, upper = upper)
    } else {
      out <- scale * crps_cnorm(y / scale, lower = lower, upper = upper)
      ind <- scale == 0 & lower <= upper
      out[ind] <-
        rep_len(abs(y - pmax(lower, 0) - pmin(upper, 0)), length(out))[ind]
      out
    }
  }
}


# truncated
#' @rdname scores_norm
#' @usage NULL
#' @export
crps_tnorm <- function(y, location = 0, scale = 1,
                       lower = -Inf, upper = Inf) {
  if (!identical(location, 0)) {
    y <- y - location
    if (!identical(lower, -Inf)) lower <- lower - location
    if (!identical(upper,  Inf)) upper <- upper - location
  }
  
  if (identical(scale, 1)) {
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
  } else {
    scale[scale < 0] <- NaN
    if (!identical(lower, -Inf)) lower <- lower / scale
    if (!identical(upper,  Inf)) upper <- upper / scale
    if (all(scale > 0, na.rm = TRUE)) {
      scale * crps_tnorm(y / scale, lower = lower, upper = upper)
    } else {
      out <- scale * crps_tnorm(y / scale, lower = lower, upper = upper)
      ind <- scale == 0 & lower < 0 & upper > 0
      out[ind] <- rep_len(abs(y), length(out))[ind]
      out
    }
  }
}

# generalized truncated/censored
#' @rdname scores_norm
#' @usage NULL
#' @export
crps_gtcnorm <- function(y, location = 0, scale = 1,
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
      p_l <- pnorm(lower)
      out_l1 <- lower * lmass^2
      out_l1[lmass == 0] <- 0
      out_l2 <- 2 * dnorm(lower) * lmass
      out_l3 <- pnorm(lower, sd = sqrt(0.5))
      z <- pmax(lower, z)
    }
    if (!identical(upper, Inf) || !identical(umass, 0)) {
      umass[umass < 0 & umass > 1] <- NaN
      p_u <- pnorm(upper)
      out_u1 <- upper * umass^2
      out_u1[umass == 0] <- 0
      out_u2 <- 2 * dnorm(upper) * umass
      out_u3 <- pnorm(upper, sd = sqrt(0.5))
      z <- pmin(upper, z)
    }
    a1 <- p_u - p_l
    a2 <- 1 - (umass + lmass)
    a2[a2 < 0 | a2 > 1] <- NaN
    b <- out_u3 - out_l3
    b[b == 0] <- NaN
    
    out <- out_u1 - out_l1 +
      (z * (2 * a2 * pnorm(z) -
             (1 - 2 * lmass) * p_u -
             (1 - 2 * umass) * p_l) +
         (2 * dnorm(z) - out_u2 - out_l2 -
            a2 * b / a1 / sqrt(pi)
         ) * a2
      ) / a1

    out[lower > upper] <- NaN
    out[lower == upper] <- 0
    out + abs(y - z)
  } else {
    scale[scale < 0] <- NaN
    if (!identical(lower, -Inf)) lower <- lower / scale
    if (!identical(upper,  Inf)) upper <- upper / scale
    if (all(scale > 0, na.rm = TRUE)) {
      scale * crps_gtcnorm(y / scale,
                           lower = lower, upper = upper,
                           lmass = lmass, umass = umass)
    } else {
      out <- scale * crps_gtcnorm(y / scale,
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

### log score ###

#' @rdname scores_norm
#' @usage NULL
#' @export
logs_norm <- function(y, mean = 0, sd = 1, location = mean, scale = sd) {
  if (!missing(mean) && !missing(location))
    stop("specify 'mean' or 'location' but not both")
  if (!missing(sd) && !missing(scale))
    stop("specify 'sd' or 'scale' but not both")
  -dnorm(y, location, scale, log = TRUE)
}

#' @rdname scores_norm
#' @usage NULL
#' @export
logs_tnorm <- function(y, location = 0, scale = 1,
                    lower = -Inf, upper = Inf) {
  -fnorm(y, location, scale, lower, upper, log = TRUE)
}


### gradient (location, scale) ###

# standard
#' @rdname scores_norm
#' @usage NULL
#' @export
gradcrps_norm <- function(y, location = 0, scale = 1) {
  if (identical(location, 0) &&
      identical(scale, 1)) {
    
    term0 <- crps_norm(y)
    term1 <- 1 - 2 * pnorm(y)
    
    cbind(dloc = term1, dscale = term0 + y * term1)
  } else if (all(is.finite(scale) & scale > 0)) {
    gradcrps_norm((y - location) / scale)
  } else {
    input <- data.frame(z = y - location, scale = scale)
    out <- matrix(NaN, dim(input)[1L], 2,
                  dimnames = list(NULL, c("dloc", "dscale")))
    isNaN <- with(input, is.na(scale) | scale <= 0)
    ind2 <- !isNaN
    if (any(ind2)) {
      out[ind2] <- with(input[ind2, ],
                        gradcrps_norm(z / scale))
    }
    out
  }
}


# censored
#' @rdname scores_norm
#' @usage NULL
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
#' @rdname scores_norm
#' @usage NULL
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
#' @rdname scores_norm
#' @usage NULL
#' @export
hesscrps_norm <- function(y , location = 0, scale = 1) {
  if (identical(location, 0) &&
      identical(scale, 1)) {
    
    term1 <- dnorm(y)
    
    d2loc <- term1
    dloc.dscale <- dscale.dloc <- term1 * y
    d2scale <- term1 * y^2
    
    2 * cbind(d2loc, d2scale, dloc.dscale, dscale.dloc)
  } else if (all(is.finite(scale) & scale > 0)) {
    hesscrps_norm((y - location) / scale) / scale
  } else {
    input <- data.frame(z = y - location, scale = scale)
    out <- matrix(NaN, dim(input)[1L], 4,
                  dimnames = list(NULL, c("d2loc", "d2scale",
                                          "dloc.dscale", "dscale.dloc")))
    isNaN <- with(input, is.na(scale) | scale <= 0)
    ind2 <- !isNaN
    if (any(ind2)) {
      out[ind2] <- with(input[ind2, ],
                        hesscrps_norm(z / scale) / scale)
    }
    out
  }
}


# censored
#' @rdname scores_norm
#' @usage NULL
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
#' @rdname scores_norm
#' @usage NULL
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

################################## Checks ######################################
check_crps_norm <- function(input) {
  required <- list(c("y", "location", "scale"),
                   c("y", "mean", "sd"))
  checkNames2(required, names(input))
  checkNumeric(input)
  checkVector(input)
  
  if ("scale" %in% names(input)) {
    if (any(input$scale <= 0))
      stop("Parameter 'scale' contains non-positive values.")
  }
  if ("sd" %in% names(input)) {
    if (any(input$sd <= 0))
      stop("Parameter 'sd' contains non-positive values.")
  }
}

check_crps_cnorm <- function(input) {
  required <- c("y", "location", "scale", "lower", "upper")
  checkNames1(required, names(input))
  checkNumeric(input, infinite_exception = c("lower", "upper"))
  checkVector(input)
  
  if (any(input$scale <= 0))
    stop("Parameter 'scale' contains non-positive values.")
  if (any(input$lower > input$upper))
    stop("Parameter 'lower' contains values greater than corresponding values in 'upper'.")
}

check_crps_tnorm <- check_crps_cnorm

check_crps_gtcnorm <- function(input) {
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

check_logs_norm <- check_crps_norm

check_logs_tnorm <- check_crps_tnorm
