#' Calculating scores for Student's \eqn{t}-distribution
#' 
#' These functions calculate scores (CRPS, logarithmic score) and their gradient and Hessian with respect
#' to the parameters of a location-scale transformed Student's
#' \eqn{t}-distribution. Furthermore, the censoring transformation and
#' the truncation transformation may be introduced on top of the
#' location-scale transformed normal distribution.
#' 
#' @usage
#' ## score functions
#' crps_t(y, df, location = 0, scale = 1)
#' crps_ct(y, df, location = 0, scale = 1, lower = -Inf, upper = Inf)
#' crps_tt(y, df, location = 0, scale = 1, lower = -Inf, upper = Inf)
#' crps_gtct(y, df, location = 0, scale = 1, lower = -Inf, upper = Inf, lmass = 0, umass = 0)
#' logs_t(y, df, location = 0, scale = 1)
#' logs_tt(y, df, location = 0, scale = 1, lower = -Inf, upper = Inf)
#' dss_t(y, df, location = 0, scale = 1)
#'
#' ## gradient (location, scale) functions
#' gradcrps_t(y, df, location = 0, scale = 1)
#' gradcrps_ct(y, df, location = 0, scale = 1, lower = -Inf, upper = Inf)
#' gradcrps_tt(y, df, location = 0, scale = 1, lower = -Inf, upper = Inf)
#'
#' ## Hessian (location, scale) functions
#' hesscrps_t(y, df, location = 0, scale = 1)
#' hesscrps_ct(y, df, location = 0, scale = 1, lower = -Inf, upper = Inf)
#' hesscrps_tt(y, df, location = 0, scale = 1, lower = -Inf, upper = Inf)
#' 
#' @param y vector of observations.
#' @param df vector of degrees of freedom.
#' @param location vector of location parameters.
#' @param scale vector of scale paramters.
#' @param lower,upper lower and upper truncation/censoring bounds.
#' @param lmass,umass vectors of point masses in \code{lower} and \code{upper}
#'  respectively. 
#' @return For the CRPS functions: a vector of score values.
#' 
#' For the gradient and Hessian functions: a matrix with column names
#' corresponding to the respective partial derivatives.
#' @name scores_t
NULL


### crps ###

# standard
#' @rdname scores_t
#' @usage NULL
#' @export
crps_t <- function(y, df, location = 0, scale = 1) {
  if (!identical(location, 0)) y <- y - location
  
  if (identical(scale, 1)) {
    df[df <= 1] <- NaN
    bfrac <- beta(0.5, df - 0.5) / beta(0.5, 0.5 * df)^2
    y * (2 * pt(y, df) - 1) + 2 / (df - 1) *
      (dt(y, df) * (df + y^2) - sqrt(df) * bfrac)
  } else {
    scale[scale < 0] <- NaN
    if (all(scale > 0, na.rm = TRUE)) {
      scale * crps_t(y / scale, df)
    } else {
      out <- scale * crps_t(y / scale, df)
      ind1 <- df == Inf
      ind2 <- scale == 0
      out[ind1] <- rep_len(scale * crps_norm(y / scale), length(out))[ind1]
      out[ind2] <- rep_len(abs(y), length(out))[ind2]
      out
    }
  }
}


# censored
#' @rdname scores_t
#' @usage NULL
#' @export
crps_ct <- function(y, df, location = 0, scale = 1,
                    lower = -Inf, upper = Inf) {
  if (!identical(location, 0)) {
    y <- y - location
    if (!identical(lower, -Inf)) lower <- lower - location
    if (!identical(upper,  Inf)) upper <- upper - location
  }
  
  if (identical(scale, 1)) {
    df[df <= 1] <- NaN
    
    out_l1 <- out_l2 <- out_u1 <- 0
    out_u2 <- 1
    z <- y
    if (!identical(lower, -Inf)) {
      p_l <- pt(lower, df)
      G_l <- -(df + lower^2) / (df - 1) * dt(lower, df)
      pb_l <- pbeta(df / (df + lower^2), df - 0.5, 0.5)
      out_l1 <- -lower * p_l^2 + 2 * G_l * p_l
      out_l1[lower == -Inf] <- 0
      out_l2 <- 0.5 * ifelse(p_l <= 0.5, pb_l, 2 - pb_l)
      z <- pmax(lower, z)
    }
    if (!identical(upper, Inf)) {
      p_u <- pt(upper, df, lower.tail = FALSE)
      G_u <- -(df + upper^2) / (df - 1) * dt(upper, df)
      pb_u <- pbeta(df / (df + upper^2), df - 0.5, 0.5)
      out_u1 <- upper * p_u^2 + 2 * G_u * p_u
      out_u1[upper == Inf] <- 0
      out_u2 <- 0.5 * ifelse(p_u >= 0.5, pb_u, 2 - pb_u)
      z <- pmin(upper, z)
    }
    b <- out_u2 - out_l2
    bfrac <- 2 * sqrt(df) / (df - 1) *
      beta(0.5, df - 0.5) / beta(0.5, 0.5 * df)^2
    G_z <- -(df + z^2) / (df - 1) * dt(z, df)
    out_z <- z * (2 * pt(z, df) - 1) - 2 * G_z
    
    out <- out_z + out_l1 + out_u1 - b * bfrac
    out[lower > upper] <- NaN
    out[lower == upper] <- 0
    out + abs(y - z)
  } else {
    scale[scale < 0] <- NaN
    if (!identical(lower, -Inf)) lower <- lower / scale
    if (!identical(upper,  Inf)) upper <- lower / scale
    if (all(scale > 0, na.rm = TRUE)) {
      scale * crps_ct(y / scale, df, lower = lower, upper = upper)
    } else {
      out <- scale * crps_ct(y / scale, df, lower = lower, upper = upper)
      ind1 <- df == Inf
      ind2 <- scale == 0 & lower <= upper
      out[ind1] <-
        rep_len(scale * crps_cnorm(y / scale, lower = lower, upper = upper),
                length(out))[ind1]
      out[ind2] <-
        rep_len(abs(y - pmax(lower, 0) - pmin(upper, 0)), length(out))[ind2]
      out
    }
  }
}


# truncated
#' @rdname scores_t
#' @usage NULL
#' @export
crps_tt <- function(y, df, location = 0, scale = 1,
                    lower = -Inf, upper = Inf) {
  if (!identical(location, 0)) {
    y <- y - location
    if (!identical(lower, -Inf)) lower <- lower - location
    if (!identical(upper,  Inf)) upper <- upper - location
  }
  
  if (identical(scale, 1)) {
    df[df <= 1] <- NaN
    
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
    }
    
    out_l <- p_l <- 0
    out_u <- p_u <- 1
    z <- y
    if (!identical(lower, -Inf)) {
      p_l <- pt(lower, df)
      pb_l <- pbeta(df / (df + lower^2), df - 0.5, 0.5)
      out_l <- 0.5 * ifelse(p_l <= 0.5, pb_l, 2 - pb_l)
      z <- pmax(lower, z)
    }
    if (!identical(upper, Inf)) {
      p_u <- pt(upper, df)
      pb_u <- pbeta(df / (df + upper^2), df - 0.5, 0.5)
      out_u <- 0.5 * ifelse(p_u <= 0.5, pb_u, 2 - pb_u)
      z <- pmin(upper, z)
    }
    a <- p_u - p_l
    b <- out_u - out_l
    b[b == 0] <- NaN
    bfrac <- 2 * sqrt(df) / (df - 1) *
      beta(0.5, df - 0.5) / beta(0.5, 0.5 * df)^2
    G_z <- -(df + z^2) / (df - 1) * dt(z, df) 
    out_z <- z * (2 * pt(z, df) - p_l - p_u) - 2 * G_z
    
    out <- (out_z - b / a * bfrac) / a
    out[lower > upper] <- NaN
    out[lower == upper] <- 0
    out + abs(y - z)
  } else {
    scale[scale < 0] <- NaN
    if (!identical(lower, -Inf)) lower <- lower / scale
    if (!identical(upper,  Inf)) upper <- lower / scale
    if (all(scale > 0, na.rm = TRUE)) {
      scale * crps_tt(y / scale, df, lower = lower, upper = upper)
    } else {
      out <- scale * crps_tt(y / scale, df, lower = lower, upper = upper)
      ind1 <- df == Inf
      ind2 <- scale == 0 & lower < 0 & upper > 0
      out[ind1] <-
        rep_len(scale * crps_tnorm(y / scale, lower = lower, upper = upper),
                length(out))[ind1]
      out[ind2] <- rep_len(abs(y), length(out))[ind2]
      out
    }
  }
}


# generalized truncated/censored
#' @rdname scores_t
#' @usage NULL
#' @export
crps_gtct <- function(y, df, location = 0, scale = 1,
                      lower = -Inf, upper = Inf,
                      lmass = 0, umass = 0) {
  if (!identical(location, 0)) {
    y <- y - location
    if (!identical(lower, -Inf)) lower <- lower - location
    if (!identical(upper,  Inf)) upper <- upper - location
  }
  
  if (identical(scale, 1)) {
    df[df <= 1] <- NaN
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
      p_l <- pt(lower, df)
      G_l <- -(df + lower^2) / (df - 1) * dt(lower, df)
      G_l[lower == -Inf] <- 0
      pb_l <- pbeta(df / (df + lower^2), df - 0.5, 0.5)
      out_l1 <- lower * lmass^2
      out_l1[lmass == 0] <- 0
      out_l2 <- 2 * G_l * lmass
      out_l3 <- 0.5 * ifelse(p_l <= 0.5, pb_l, 2 - pb_l)
      z <- pmax(lower, z)
    }
    if (!identical(upper, Inf) || !identical(umass, 0)) {
      umass[umass < 0 & umass > 1] <- NaN
      p_u <- pt(upper, df)
      G_u <- -(df + upper^2) / (df - 1) * dt(upper, df)
      G_u[upper == Inf] <- 0
      pb_u <- pbeta(df / (df + upper^2), df - 0.5, 0.5)
      out_u1 <- upper * umass^2
      out_u1[umass == 0] <- 0
      out_u2 <- 2 * G_u * umass
      out_u3 <- 0.5 * ifelse(p_u <= 0.5, pb_u, 2 - pb_u)
      z <- pmin(upper, z)
    }
    a1 <- p_u - p_l
    a2 <- 1 - (umass + lmass)
    a2[a2 < 0 | a2 > 1] <- NaN
    b <- out_u3 - out_l3
    b[b == 0] <- NaN
    bfrac <- 2 * sqrt(df) / (df - 1) *
      beta(0.5, df - 0.5) / beta(0.5, 0.5 * df)^2
    G_z <- -(df + z^2) / (df - 1) * dt(z, df) 
    
    out <- out_u1 - out_l1 +
      (z * (2 * a2 * pt(z, df) -
              (1 - 2 * lmass) * p_u -
              (1 - 2 * umass) * p_l) -
         (2 * G_z - out_u2 - out_l2 +
            a2 * b / a1 * bfrac
         ) * a2
      ) / a1
    
    out[lower > upper] <- NaN
    out[lower == upper] <- 0
    out + abs(y - z)
  } else {
    scale[scale < 0] <- NaN
    if (!identical(lower,-Inf)) lower <- lower / scale
    if (!identical(upper,  Inf)) upper <- upper / scale
    if (all(scale > 0, na.rm = TRUE)) {
      scale * crps_gtct(y / scale, df,
                        lower = lower, upper = upper,
                        lmass = lmass, umass = umass)
    } else {
      out <- scale * crps_gtct(y / scale, df,
                               lower = lower, upper = upper,
                               lmass = lmass, umass = umass)
      ind1 <- df == Inf
      ind2 <- scale == 0 & lower < 0 & upper > 0
      out[ind1] <-
        rep_len(scale * crps_gtcnorm(y / scale,
                                     lower = lower, upper = upper,
                                     lmass = lmass, umass = umass),
                length(out))[ind1]
      out[ind2] <-
        rep_len(
          (pmin(y, 0) - lower) * lmass ^ 2 - pmin(y, 0) * (1 - lmass) ^ 2 +
            (upper - pmax(y, 0)) * umass ^ 2 + pmax(y, 0) * (1 - umass) ^ 2,
          length(out)
        )[ind2]
      out
    }
  }
}

#' @rdname scores_t
#' @usage NULL
#' @export
logs_t <- function(y, df, location = 0, scale = 1) {
  -ft(y, df, location, scale, log = TRUE)
}

#' @rdname scores_t
#' @usage NULL
#' @export
logs_tt <- function(y, df, location = 0, scale = 1,
                    lower = -Inf, upper = Inf) {
  -ft(y, df, location, scale, lower, upper, log = TRUE)
}


#' @rdname scores_t
#' @usage NULL
#' @export
dss_t <- function(y, df, location = 0, scale = 1) {
  if (!identical(location, 0)) y <- y - location
  df[df <= 2] <- NaN
  scale[scale <= 0] <- NaN
  v <- scale^2 * ifelse(is.infinite(df), 1, df / (df - 2))
  y^2 / v  + log(v)
}


### gradient (location, scale) ###

# standard
#' @rdname scores_t
#' @usage NULL
#' @export
gradcrps_t <- function(y , df, location = 0, scale = 1) {
  all_df_in_1_to_Inf <- all(is.finite(df) & df > 1)
  if (all_df_in_1_to_Inf &&
      identical(location, 0) &&
      identical(scale, 1)) {
    
    term0 <- crps_t(y, df)
    term1 <- 1 - 2 * pt(y, df)
    
    cbind(dloc = term1, dscale = term0 + y * term1)
  } else if (all_df_in_1_to_Inf &&
             all(is.finite(scale) & scale > 0)) {
    gradcrps_t((y - location) / scale, df)
  } else {
    input <- data.frame(z = y - location,
                        df = df,
                        scale = scale)
    out <- rep(NaN, dim(input)[1L])
    isNaN <- with(input, is.na(df) | df <= 1 |
                    is.na(scale) | scale <= 0)
    ind2 <- !isNaN & input$df == Inf
    ind3 <- !isNaN & !ind2
    if (any(ind2)) {
      out[ind2] <- with(input[ind2, ], gradcrps_norm(z / scale))
    }
    if (any(ind3)) {
      out[ind3] <- with(input[ind3, ], gradcrps_t(z / scale, df))
    }
    out
  }
}

# censored
#' @rdname scores_t
#' @usage NULL
#' @export
gradcrps_ct <- function(y, df, location = 0, scale = 1,
                        lower = -Inf, upper = Inf) {
  nan_in_bounds <- anyNA(lower) || anyNA(upper)
  all_df_in_1_to_Inf <- all(is.finite(df) & df > 1)
  if (!nan_in_bounds &&
      all_df_in_1_to_Inf &&
      identical(location, 0) &&
      identical(scale, 1)) {
    z <- y
    if (!identical(lower, -Inf)) z <- pmax(lower, z)
    if (!identical(upper,  Inf)) z <- pmin(upper, z)
    
    term0 <- crps_ct(z, df, lower = lower, upper = upper)
    term1 <- 2 * pt(z, df) - 1
    term2 <- pt(lower, df)^2
    term3 <- pt(upper, df, lower.tail = FALSE)^2
    
    dloc <- -term1 + term2 - term3
    dscale <- term0 - term1 * z +
      term2 * ifelse(lower == -Inf, 0, lower) -
      term3 * ifelse(upper ==  Inf, 0, upper)
    
    out <- cbind(dloc, dscale)
    out[lower > upper, ] <- NaN
    out[lower == upper, ] <- 0
    out
  } else if (!nan_in_bounds &&
             all_df_in_1_to_Inf &&
             all(is.finite(scale) & scale > 0)) {
    if (!identical(lower, -Inf)) lower <- (lower - location) / scale
    if (!identical(upper,  Inf)) upper <- (upper - location) / scale
    gradcrps_ct((y - location) / scale, df,
                lower = lower, upper = upper)
  } else {
    input <- data.frame(z = y - location, scale = scale,
                        lower = lower - location,
                        upper = upper - location)
    out <- matrix(NaN, dim(input)[1L], 2,
                  dimnames = list(NULL, c("dloc", "dscale")))
    isNaN <- with(input, {
      is.na(df) | df <= 1 |
        is.na(scale) | scale <= 0 |
        is.na(lower) | is.na(upper)
    })
    ind2 <- !isNaN & is.infinite(df)
    ind3 <- !isNaN & !ind2
    if (any(ind2)) {
      out[ind2, ] <- with(input[ind2, ], {
        gradcrps_cnorm(z / scale,
                       lower = lower / scale,
                       upper = upper / scale)
      })
    }
    if (any(ind3)) {
      out[ind3, ] <- with(input[ind3, ], {
        gradcrps_ct(z / scale, df,
                    lower = lower / scale,
                    upper = upper / scale)
      })
    }
    out
  }
}


# truncated
#' @rdname scores_t
#' @usage NULL
#' @export
gradcrps_tt <- function(y, df, location = 0, scale = 1,
                        lower = -Inf, upper = Inf) {
  nan_in_bounds <- anyNA(lower) || anyNA(upper)
  all_df_in_1_to_Inf <- all(is.finite(df) & df > 1)
  
  if (!nan_in_bounds &&
      all_df_in_1_to_Inf &&
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
      p_l <- pt(lower, df)
      d_l <- dt(lower, df)
      G_l <- (df + lower^2) / (df - 1) * d_l
      G_l[is.infinite(lower)] <- 0
      z <- pmax(lower, z)
    }
    if (!identical(upper, Inf)) {
      p_u <- pt(upper, df)
      d_u <- dt(upper, df)
      G_u <- (df + upper^2) / (df - 1) * d_u
      G_u[is.infinite(upper)] <- 0
      z <- pmin(upper, z)
    }
    a <- p_u - p_l
    p_z <- pt(z, df)
    d_z <- dt(z, df)
    G_z <- (df + z^2) / (df - 1) * d_z
    
    term0 <- crps_tt(z, df, lower = lower, upper = upper)
    term1 <- (p_l + p_u - 2 * p_z) / a
    term2 <- 2 * d_l / a * 
      ((z * p_l - z * p_z - G_l + G_z) / a + term0)
    term3 <- 2 * d_u / a *
      ((z * p_u - z * p_z - G_u + G_z) / a + term0)
    
    dloc <- (term1 - term2 + term3) * sign
    dscale <- term0 + term1 * z - 
      term2 * ifelse(lower == -Inf, 0, lower) +
      term3 * ifelse(upper ==  Inf, 0, upper)
    
    out <- cbind(dloc, dscale)
    out[lower > upper, ] <- NaN
    out[lower == upper, ] <- 0
    out
  } else if (!nan_in_bounds &&
             all_df_in_1_to_Inf &&
             all(is.finite(scale) & scale > 0)) {
    if (!identical(lower, -Inf)) lower <- (lower - location) / scale
    if (!identical(upper,  Inf)) upper <- (upper - location) / scale
    gradcrps_tt((y - location) / scale, df,
                lower = lower, upper = upper)
  } else {
    input <- data.frame(z = y - location, df = df, scale = scale,
                        lower = lower - location,
                        upper = upper - location)
    out <- matrix(NaN, dim(input)[1L], 2,
                  dimnames = list(NULL, c("dloc", "dscale")))
    isNaN <- with(input, {
      is.na(df) | df <= 1 |
        is.na(scale) | scale <= 0 |
        is.na(lower) | is.na(upper)
    })
    ind2 <- !isNaN & is.infinite(df)
    ind3 <- !isNaN & !ind2
    if (any(ind2)) {
      out[ind2, ] <- with(input[ind2, ], {
        gradcrps_tnorm(z / scale,
                       lower = lower / scale,
                       upper = upper / scale)
      })
    }
    if (any(ind3)) {
      out[ind3, ] <- with(input[ind3, ], {
        gradcrps_tt(z / scale, df,
                    lower = lower / scale,
                    upper = upper / scale)
      })
    }
    out
  }
}

### Hessian (location, scale) ###

# standard
#' @rdname scores_t
#' @usage NULL
#' @export
hesscrps_t <- function(y , df, location = 0, scale = 1) {
  all_df_in_1_to_Inf <- all(is.finite(df) & df > 1)
  if (all_df_in_1_to_Inf &&
      identical(location, 0) &&
      identical(scale, 1)) {
    
    term1 <- dt(y, df)
    d2loc <- term1
    dloc.dscale <- dscale.dloc <- term1 * y
    d2scale <- term1 * y^2
    
    2 * cbind(d2loc, d2scale, dloc.dscale, dscale.dloc)
  } else if (all_df_in_1_to_Inf &&
             all(is.finite(scale) & scale > 0)) {
    hesscrps_t((y - location) / scale, df) / scale
  } else {
    input <- data.frame(z = y - location,
                        df = df,
                        scale = scale)
    out <- rep(NaN, dim(input)[1L], 4,
               dimnames = list(NULL, c("d2loc", "d2scale",
                                       "dloc.dscale", "dscale.dloc")))
    isNaN <- with(input, is.na(df) | df <= 1 |
                    is.na(scale) | scale <= 0)
    ind2 <- !isNaN & input$df == Inf
    ind3 <- !isNaN & !ind2
    if (any(ind2)) {
      out[ind2] <- with(input[ind2, ], hesscrps_norm(z / scale) / scale)
    }
    if (any(ind3)) {
      out[ind3] <- with(input[ind3, ], gesscrps_t(z / scale, df) / scale)
    }
    out
  }
}


# censored
#' @rdname scores_t
#' @usage NULL
#' @export
hesscrps_ct <- function(y, df, location = 0, scale = 1,
                        lower = -Inf, upper = Inf) {
  nan_in_bounds <- anyNA(lower) || anyNA(upper)
  all_df_in_1_to_Inf <- all(is.finite(df) & df > 1)
  if (!nan_in_bounds &&
      all_df_in_1_to_Inf &&
      identical(location, 0) &&
      identical(scale, 1)) {
    z <- y
    if (!identical(lower, -Inf)) z <- pmax(lower, z)
    if (!identical(upper,  Inf)) z <- pmin(upper, z)
    
    term1 <- dt(z, df)
    term2 <- dt(lower, df) * pt(lower, df)
    term3 <- dt(upper, df) * pt(upper, df, lower.tail = FALSE)
    
    d2mu <- term1 - term2 - term3
    dmu.dsigma <- dsigma.dmu <- term1 * z -
      term2 * ifelse(is.finite(lower), lower, 0) -
      term3 * ifelse(is.finite(upper), upper, 0)
    d2sigma <- term1 * z^2 -
      term2 * ifelse(is.finite(lower), lower^2, 0) -
      term3 * ifelse(is.finite(upper), upper^2, 0)
    
    out <- 2 * cbind(d2mu, d2sigma, dmu.dsigma, dsigma.dmu)
    out[lower > upper, ] <- NaN
    out[lower == upper, ] <- 0
    out
  } else if (!nan_in_bounds &&
             all_df_in_1_to_Inf &&
             all(is.finite(scale) & scale > 0)) {
    if (!identical(lower, -Inf)) lower <- (lower - location) / scale
    if (!identical(upper,  Inf)) upper <- (upper - location) / scale
    hesscrps_ct((y - location) / scale, df,
                lower = lower, upper = upper)
  } else {
    input <- data.frame(z = y - location, scale = scale,
                        lower = lower - location,
                        upper = upper - location)
    out <- rep(NaN, dim(input)[1L], 4,
               dimnames = list(NULL, c("d2loc", "d2scale",
                                       "dloc.dscale", "dscale.dloc")))
    isNaN <- is.na(input$z) |
      is.na(input$df) | input$df <= 1 |
      is.na(input$scale) | input$scale <= 0
    ind2 <- input$df == Inf & !isNaN
    ind3 <- !isNaN & !ind2
    if (any(ind2)) {
      out[ind2] <- with(input[ind2, ],
                        scale * hesscrps_cnorm(z / scale,
                                               lower = lower / scale,
                                               upper = upper / scale))
    }
    if (any(ind3)) {
      out[ind3] <- with(input[ind3, ],
                        scale * hesscrps_ct(z / scale, df,
                                            lower = lower / scale,
                                            upper = upper / scale))
    }
    out
  }
}


# truncated
#' @rdname scores_t
#' @usage NULL
#' @export
hesscrps_tt <- function(y, df, location = 0, scale = 1,
                        lower = -Inf, upper = Inf) {
  nan_in_bounds <- anyNA(lower) || anyNA(upper)
  all_df_in_1_to_Inf <- all(is.finite(df) & df > 1)
  if (!nan_in_bounds &&
      all_df_in_1_to_Inf &&
      identical(location, 0) &&
      identical(scale, 1)) {
    
    z <- y
    if (!identical(lower, -Inf)) z <- pmax(lower, z)
    if (!identical(upper,  Inf)) z <- pmin(upper, z)
    
    z_lb_ub <- cbind(z, lb = lower, ub = upper)
    
    CDF <- pt(z_lb_ub, df)
    denom <- CDF[, "ub"] - CDF[, "lb"]
    CDF <- CDF / denom
    PDF <- dt(z_lb_ub, df) / denom
    PDFp_over_PDF <- -(df + 1) * z_lb_ub / (df + z_lb_ub^2)
    PDFp_over_PDF[is.infinite(z_lb_ub)] <- 0
    G <- -(df + z_lb_ub^2) / (df - 1) * PDF
    G[is.infinite(z_lb_ub)] <- 0
    CRPS <- crps_tt(z, df, lower = lower, upper = upper)
    
    out <- calcHess_trunc(z, 1, lower, upper,
                          CDF, PDF, PDFp_over_PDF, G, CRPS)
    out[lower > upper, ] <- NaN
    out[lower == upper, ] <- 0
    out
  } else if (!nan_in_bounds &&
             all_df_in_1_to_Inf &&
             all(is.finite(scale) & scale > 0)) {
    if (!identical(lower, -Inf)) lower <- (lower - location) / scale
    if (!identical(upper,  Inf)) upper <- (upper - location) / scale
    hesscrps_tt((y - location) / scale, df,
                lower = lower, upper = upper)
  } else {
    input <- data.frame(z = y - location, scale = scale,
                        lower = lower - location,
                        upper = upper - location)
    out <- rep(NaN, dim(input)[1L], 4,
               dimnames = list(NULL, c("d2loc", "d2scale",
                                       "dloc.dscale", "dscale.dloc")))
    isNaN <- is.na(input$z) |
      is.na(input$df) | input$df <= 1 |
      is.na(input$scale) | input$scale <= 0
    ind2 <- input$df == Inf & !isNaN
    ind3 <- !isNaN & !ind2
    if (any(ind2)) {
      out[ind2] <- with(input[ind2, ],
                        scale * hesscrps_tnorm(z / scale,
                                               lower = lower / scale,
                                               upper = upper / scale))
    }
    if (any(ind3)) {
      out[ind3] <- with(input[ind3, ],
                        scale * hesscrps_tt(z / scale, df,
                                            lower = lower / scale,
                                            upper = upper / scale))
    }
    out
  }
}

################################## Checks ######################################
check_crps_t <- function(input) {
  required <- c("y", "df", "location", "scale")
  checkNames1(required, names(input))
  checkNumeric(input, infinite_exception = "df")
  checkVector(input)
  
  if (any(input$scale <= 0))
    stop("Parameter 'scale' contains non-positive values.")
  if (any(input$df <= 1)) {
    stop(paste("Parameter 'df' contains values less than or equal to 1.",
               "The CRPS does not exist."))
  }
}

check_crps_ct <- function(input) {
  required <- c("y", "location", "scale", "lower", "upper")
  checkNames1(required, names(input))
  checkNumeric(input, infinite_exception = c("lower", "upper"))
  checkVector(input)
  
  if (any(input$scale <= 0))
    stop("Parameter 'scale' contains non-positive values.")
  if (any(input$df <= 1)) {
    stop(paste("Parameter 'df' contains values less than or equal to 1.",
               "The CRPS does not exist."))
  }
  if (any(input$lower > input$upper))
    stop("Parameter 'lower' contains values greater than corresponding values in 'upper'.")
}

check_crps_tt <- check_crps_ct

check_crps_gtclogis <- function(input) {
  required <- c("y", "location", "scale", "lower", "upper", "lmass", "umass")
  checkNames1(required, names(input))
  checkNumeric(input, infinite_exception = c("lower", "upper"))
  checkVector(input)
  
  if (any(input$scale <= 0))
    stop("Parameter 'scale' contains non-positive values.")
  if (any(input$df <= 1)) {
    stop(paste("Parameter 'df' contains values less than or equal to 1.",
               "The CRPS does not exist."))
  }
  if (any(input$lower > input$upper))
    stop("Parameter 'lower' contains values greater than corresponding values in 'upper'.")
  if (any(input$lmass < 0 | input$lmass > 1))
    stop("Parameter 'lmass' contains values not in [0, 1].")
  if (any(input$umass < 0 | input$umass > 1))
    stop("Parameter 'umass' contains values not in [0, 1].")
  if (any(input$umass + input$lmass > 1))
    stop("Values in 'lmass' and 'umass' add up to more than 1.")
}

check_logs_t <- function(input) {
  required <- c("y", "df", "location", "scale")
  checkNames1(required, names(input))
  checkNumeric(input, infinite_exception = "df")
  checkVector(input)
  
  if (any(input$scale <= 0))
    stop("Parameter 'scale' contains non-positive values.")
  if (any(input$df <= 0))
    stop("Parameter 'df' contains non-positive values.")
}

check_logs_tt <- function(input) {
  required <- c("y", "location", "scale", "lower", "upper")
  checkNames1(required, names(input))
  checkNumeric(input, infinite_exception = c("lower", "upper"))
  checkVector(input)
  
  if (any(input$scale <= 0))
    stop("Parameter 'scale' contains non-positive values.")
  if (any(input$df <= 0))
    stop("Parameter 'df' contains non-positive values.")
  if (any(input$lower > input$upper))
    stop("Parameter 'lower' contains values greater than corresponding values in 'upper'.")
}
