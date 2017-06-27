### crps ###

# standard
#' @export
crps_t <- function(y, df, location = 0, scale = 1) {
  all_df_in_1_to_Inf <- all(is.finite(df) & df > 1)
  if (all_df_in_1_to_Inf &&
      identical(location, 0) &&
      identical(scale, 1)) {
    bfrac <- beta(0.5, df - 0.5) / beta(0.5, 0.5 * df)^2
    y * (2 * pt(y, df) - 1) + 2 / (df - 1) *
      (dt(y, df) * (df + y^2) - sqrt(df) * bfrac)
  } else if (all_df_in_1_to_Inf &&
             all(is.finite(scale) & scale > 0)) {
    scale * crps_t((y - location) / scale, df)
  } else {
    input <- data.frame(z = abs(y - location),
                        df = df,
                        scale = scale)
    out <- rep(NaN, dim(input)[1L])
    isNaN <- with(input, is.na(df) | df <= 1 |
                    is.na(scale) | scale < 0)
    ind1 <- input$scale == 0 & !isNaN
    ind2 <- input$df == Inf & !isNaN & !ind1
    ind3 <- !isNaN & !ind1 & !ind2
    if (any(ind1)) {
      out[ind1] <- with(input[ind1, ], z)
    }
    if (any(ind2)) {
      out[ind2] <- with(input[ind2, ], crps_norm(-z / scale))
    }
    if (any(ind3)) {
      out[ind3] <- with(input[ind3, ], crps_t(-z / scale, df))
    }
    out
  }
}


# censored
#' @export
crps_ct <- function(y, df, location = 0, scale = 1,
                    lower = -Inf, upper = Inf) {
  nan_in_bounds <- anyNA(lower) || anyNA(upper)
  all_df_in_1_to_Inf <- all(is.finite(df) & df > 1)
  
  if (!nan_in_bounds &&
      all_df_in_1_to_Inf &&
      identical(location, 0) &&
      identical(scale, 1)) {
    
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
  } else if (!nan_in_bounds &&
             all_df_in_1_to_Inf &&
             all(is.finite(scale) & scale > 0)) {
    if (!identical(lower, -Inf)) lower <- (lower - location) / scale
    if (!identical(upper,  Inf)) upper <- (upper - location) / scale
    scale * crps_ct((y - location) / scale, df,
                    lower = lower, upper = upper)
  } else {
    input <- data.frame(z = y - location, df = df, scale = scale,
                        lower = lower - location,
                        upper = upper - location)
    
    out <- rep(NaN, dim(input)[1L])
    isNaN <- with(input, {
      is.na(input$df) | input$df <= 1 |
        is.na(input$scale) | input$scale < 0 |
        is.na(lower) | is.na(upper)
    })
    ind1 <- !isNaN & with(input, scale == 0 & lower <= upper)
    ind2 <- !isNaN & !ind1 & is.infinite(input$df)
    ind3 <- !isNaN & !ind2 & scale > 0
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
    if (any(ind3)) {
      out[ind3] <- with(input[ind3, ], {
        scale * crps_t(z / scale, df,
                       lower = lower / scale,
                       upper = upper / scale)
      })
    }
    out
  }
}


# truncated
#' @export
crps_tt <- function(y, df, location = 0, scale = 1,
                    lower = -Inf, upper = Inf) {
  nan_in_bounds <- anyNA(lower) || anyNA(upper)
  all_df_in_1_to_Inf <- all(is.finite(df) & df > 1)
  
  if (!nan_in_bounds &&
      all_df_in_1_to_Inf &&
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
  } else if (!nan_in_bounds &&
             all_df_in_1_to_Inf &&
             all(is.finite(scale) & scale > 0)) {
    if (!identical(lower, -Inf)) lower <- (lower - location) / scale
    if (!identical(upper,  Inf)) upper <- (upper - location) / scale
    scale * crps_tt((y - location) / scale, df,
                    lower = lower, upper = upper)
  } else {
    input <- data.frame(z = y - location, df = df, scale = scale,
                        lower = lower - location,
                        upper = upper - location)
    
    out <- rep(NaN, dim(input)[1L])
    isNaN <- with(input, {
      is.na(input$df) | input$df <= 1 |
        is.na(input$scale) | input$scale < 0 |
        is.na(lower) | is.na(upper)
    })
    ind1 <- !isNaN & with(input, scale == 0 & lower < 0 & upper > 0)
    ind2 <- !isNaN & !ind1 & is.infinite(input$df)
    ind3 <- !isNaN & !ind2 & scale > 0
    if (any(ind1)) out[ind1] <- abs(input$z[ind1])
    if (any(ind2)) {
      out[ind2] <- with(input[ind2, ], {
        scale * crps_tnorm(z / scale,
                           lower = lower / scale,
                           upper = upper / scale)
      })
    }
    if (any(ind3)) {
      out[ind3] <- with(input[ind3, ], {
        scale * crps_tt(z / scale, df,
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
      term2 * ifelse(lower == -Inf, lower, 0) -
      term3 * ifelse(upper ==  Inf, upper, 0)
    
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
      term2 * ifelse(lower == -Inf, lower, 0) +
      term3 * ifelse(upper ==  Inf, upper, 0)
    
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
