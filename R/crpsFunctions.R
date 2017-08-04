#### Numerical integration ####
# Note: y can be a vector, all other inputs are scalars
crps.int <- function(y, pxxx, lower, upper, rel_tol = 1e-6){
  ind <- (y > upper) - (y < lower)
  out <- numeric(length(y))
  F1 <- function(x) pxxx(x)^2
  F2 <- function(x) (1-pxxx(x))^2
  if (any(ind == -1)) {
    out[ind == -1] <- sapply(which(ind == -1), function(i) {
      s1 <- lower - y[i]
      s2 <- integrate(F2, lower, upper, rel.tol = rel_tol)$value
      s1 + s2
    })
  } else if (any(ind == 0)) {
    out[ind == 0] <- sapply(which(ind == 0), function(i) {
      s1 <- integrate(F1, lower, y[i], rel.tol = rel_tol)$value
      s2 <- integrate(F2, y[i], upper, rel.tol = rel_tol)$value
      s1 + s2
    })
  } else if (any(ind == 1)) {
    out[ind == 1] <- sapply(which(ind == 1), function(i) {
      s1 <- integrate(F1, lower, upper, rel.tol = rel_tol)$value
      s2 <- y[i] - upper
      s1 + s2
    })
  }
  return(out)
}

# mixture of normals (numerical integration)
crps.mixnorm.int <- function(y, m, s, w, rel_tol){
  Fmix <- function(z){
    sapply(z, function(r) sum(w*pnorm((r-m)/s)))
  }
  crps.int(y, Fmix, -Inf, Inf, rel_tol)
}



################################################################################
### real line

# logistic
crps.logis <- function(y, location, scale,
                       lower = -Inf, upper = Inf,
                       lmass = 0, umass = 0) {
  
  ### standard formula
  
  ind1 <- any(is.finite(lower))
  ind2 <- any(is.finite(upper))
  
  if (!ind1 & !ind2) {
    z <- y
    if (!identical(location, 0) | !identical(scale, 1)) {
      z <- (y - location) / scale
    }
    out <- z - 2 * plogis(z, log.p = TRUE) - 1
    return(scale * out)
  }
  
  ### dealing with truncation/censoring
  
  zb <- y
  if (ind1) {
    zb <- pmax(lower, zb)
    lb <- (lower - location) / scale
    
    if (is.character(lmass)) {
      n1 <- length(lb)
      n2 <- length(lmass)
      if (n1 < n2) {
        Plb <- numeric(n2)
        Plb[lmass == "cens"] <- plogis(lb)
      } else {
        Plb <- numeric(n1)
        ind <- lmass == "cens"
        Plb[ind] <- plogis(lb[ind])
      }
    } else {
      Plb <- lmass
    }
  }
  if (ind2) {
    zb <- pmin(upper, zb)
    ub <- (upper - location) / scale
    
    if (is.character(umass)) {
      n1 <- length(ub)
      n2 <- length(umass)
      if (n1 < n2) {
        Pub <- numeric(n2)
        Pub[umass == "cens"] <- plogis(ub, lower.tail = FALSE)
      } else {
        Pub <- numeric(n1)
        ind <- umass == "cens"
        Pub[ind] <- plogis(ub[ind], lower.tail = FALSE)
      }
    } else {
      Pub <- umass
    }
  }
  res <- abs(y - zb)
  zb <- (zb - location) / scale
  
  if (ind1 & ind2) {
    if (any(Plb + Pub > 1)){
      stop("Sum of 'lmass' and 'umass' exceeds 1.")
    }
    a <- ifelse(1 - Plb - Pub < 1e-12, 0,
                (1 - Plb - Pub) / (plogis(ub) - plogis(lb)))
    a[a > 1e12] <- NaN
    b <- Plb - a * plogis(lb)
    
    out_l <- a^2 * plogis(lb) +
      ifelse(a == 0, 0, a * (a + 2 * b) * plogis(-lb, log.p = TRUE)) -
      ifelse(is.finite(lb), b^2 * lb, 0)
    out_u <- a^2 * plogis(ub, lower.tail = FALSE) -
      ifelse(a == 0, 0, a * (a + 2 * b - 2) * plogis(ub, log.p = TRUE)) +
      ifelse(is.finite(ub), (a + b - 1)^2 * ub, 0)
    out_y <- zb * (2 * (a + b) - 1) - 2 * a * plogis(zb, log.p = TRUE) - a^2
  } else if (ind1 & !ind2) {
    a <- ifelse(1 - Plb < 1e-12, 0, (1 - Plb) / (1 - plogis(lb)))
    a[a > 1e12] <- NaN
    b <- Plb - a * plogis(lb)
    
    out_l <- a^2 * plogis(lb) +
      ifelse(a == 0, 0, a * (a + 2 * b) * plogis(-lb, log.p = TRUE)) -
      ifelse(is.finite(lb), b^2 * lb, 0)
    out_u <- 0
    out_y <- zb * (2 * (a + b) - 1) - 2 * a * plogis(zb, log.p = TRUE) - a^2
  } else if (!ind1 & ind2) {
    a <- ifelse(1 - Pub < 1e-12, 0, (1 - Pub) / plogis(ub))
    a[a > 1e12] <- NaN
    
    out_l <- 0
    out_u <- a^2 * plogis(ub, lower.tail = FALSE) -
      ifelse(a == 0, 0, a * (a - 2) * plogis(ub, log.p = TRUE)) +
      ifelse(is.finite(ub), (a - 1)^2 * ub, 0)
    out_y <- zb * (2 * a - 1) - 2 * a * plogis(zb, log.p = TRUE) - a^2
  }
  
  return(res + scale * (out_y + out_l + out_u))
}

# normal
crps.norm <- function(y, location, scale,
                      lower = -Inf, upper = Inf,
                      lmass = 0, umass = 0) {
  
  ### standard formula
  
  ind1 <- any(is.finite(lower))
  ind2 <- any(is.finite(upper))
  
  if (!ind1 & !ind2) {
    z <- y
    if (!identical(location, 0) | !identical(scale, 1)) {
      z <- (y - location) / scale
    }
    out <- z * (2 * pnorm(z) - 1) + 2 * dnorm(z) - 1 / sqrt(pi)
    return(scale * out)
  }
  
  ### dealing with truncation/censoring
  
  zb <- y
  if (ind1) {
    zb <- pmax(lower, zb)
    lb <- (lower - location) / scale
    
    if (is.character(lmass)) {
      n1 <- length(lb)
      n2 <- length(lmass)
      if (n1 < n2) {
        Plb <- numeric(n2)
        Plb[lmass == "cens"] <- pnorm(lb)
      } else {
        Plb <- numeric(n1)
        ind <- lmass == "cens"
        Plb[ind] <- pnorm(lb[ind])
      }
    } else {
      Plb <- lmass
    }
  }
  if (ind2) {
    zb <- pmin(upper, zb)
    ub <- (upper - location) / scale
    
    if (is.character(umass)) {
      n1 <- length(ub)
      n2 <- length(umass)
      if (n1 < n2) {
        Pub <- numeric(n2)
        Pub[umass == "cens"] <- pnorm(ub, lower.tail = FALSE)
      } else {
        Pub <- numeric(n1)
        ind <- umass == "cens"
        Pub[ind] <- pnorm(ub[ind], lower.tail = FALSE)
      }
    } else {
      Pub <- umass
    }
  }
  res <- abs(y - zb)
  zb <- (zb - location) / scale
  
  if (ind1 & ind2) {
    if (any(Plb + Pub > 1)){
      stop("Sum of 'lmass' and 'umass' exceeds 1.")
    }
    a <- ifelse(1 - Plb - Pub < 1e-12, 0,
                (1 - Plb - Pub) / (pnorm(ub) - pnorm(lb)))
    a[a > 1e12] <- NaN
    
    out_l <- -2 * a * dnorm(lb) * Plb +
      a^2 / sqrt(pi) * pnorm(lb * sqrt(2)) -
      ifelse(is.finite(lb), lb * Plb^2, 0)
    out_u <- -2 * a * dnorm(ub) * Pub +
      a^2 / sqrt(pi) * pnorm(ub * sqrt(2), lower.tail = FALSE) +
      ifelse(is.finite(ub), ub * Pub^2, 0)
    out_y <- zb * (2 * (a * (pnorm(zb) - pnorm(lb)) + Plb) - 1) +
      2 * a * dnorm(zb) - a^2 / sqrt(pi)
  } else if (ind1 & !ind2) {
    a <- ifelse(1 - Plb < 1e-12, 0, (1 - Plb) / (1 - pnorm(lb)))
    a[a > 1e12] <- NaN
    
    out_l <- -2 * a * dnorm(lb) * Plb +
      a^2 / sqrt(pi) * pnorm(lb * sqrt(2)) -
      ifelse(is.finite(lb), lb * Plb^2, 0)
    out_u <- 0
    out_y <- zb * (2 * (1 - a * pnorm(zb, lower.tail = FALSE)) - 1) +
      2 * a * dnorm(zb) - a^2 / sqrt(pi)
  } else if (!ind1 & ind2) {
    a <- ifelse(1 - Pub < 1e-12, 0, (1 - Pub) / pnorm(ub))
    a[a > 1e12] <- NaN
    
    out_l <- 0
    out_u <- -2 * a * dnorm(ub) * Pub +
      a^2 / sqrt(pi) * pnorm(ub * sqrt(2), lower.tail = FALSE) +
      ifelse(is.finite(ub), ub * Pub^2, 0)
    out_y <- zb * (2 * a * pnorm(zb) - 1) +
      2 * a * dnorm(zb) - a^2 / sqrt(pi)
  }
  
  return(res + scale * (out_y + out_l + out_u))
}

# t
crps.t <- function(y, df, location, scale,
                   lower = -Inf, upper = Inf,
                   lmass = 0, umass = 0) {
  
  if (any(!df > 1))
    stop("Parameter 'df' contains values not greater than 1. The CRPS is not defined.")
  
  ### standard formula
  
  ind <- df == Inf
  ind1 <- any(is.finite(lower))
  ind2 <- any(is.finite(upper))
  
  
  if (!ind1 & !ind2) {
    z <- y
    if (!identical(location, 0) | !identical(scale, 1)) {
      z <- (y - location) / scale
    }
    if (any(ind)) {
      if (length(z) < length(df)) z <- rep(z, len = length(df))
      
      out <- numeric(length(z))
      out[ind] <- crps.norm(z[ind], 0, 1)
      out[!ind] <- crps.t(z[!ind], df[!ind], 0, 1)
    } else {
      c1 <- z * (2 * pt(z, df) - 1)
      c2 <- dt(z, df) * (1 + z^2 / df)
      c3 <- beta(0.5, df - 0.5) / sqrt(df) / beta(0.5, 0.5 * df)^2
      out <- c1 + 2 * df / (df - 1) * (c2 - c3)
    }
    return(scale * out)
  }
  
  ### dealing with truncation/censoring
  
  zb <- y
  if (ind1) {
    zb <- pmax(lower, zb)
    lb <- (lower - location) / scale
    
    if (is.character(lmass)) {
      n <- max(length(lb), length(lmass), length(df))
      ind_cens <- lmass == "cens"
      if (any(ind_cens)) {
        Plb <- rep(pt(lb, df), len = n)
        Plb[!ind_cens] <- 0
      } else {
        Plb <- 0
      }
    } else {
      Plb <- lmass
    }
  } else {
    lb <- lower
    Plb <- lmass
  }
  if (ind2) {
    zb <- pmin(upper, zb)
    ub <- (upper - location) / scale
    
    if (is.character(umass)) {
      n <- max(length(ub), length(umass), length(df))
      ind_cens <- umass == "cens"
      if (any(ind_cens)) {
        Pub <- rep(pt(ub, df, lower.tail = FALSE), len = n)
        Pub[!ind_cens] <- 0
      } else {
        Pub <- 0
      }
    } else {
      Pub <- umass
    }
  } else {
    ub <- upper
    Pub <- umass
  }
  res <- abs(y - zb)
  zb <- (zb - location) / scale
  
  if (any(ind)) {
    if (length(zb) < length(df)) zb <- rep(zb, len = length(df))
    if (length(lb) < length(df)) lb <- rep(lb, len = length(df))
    if (length(ub) < length(df)) ub <- rep(ub, len = length(df))
    if (length(Plb) < length(df)) Plb <- rep(Plb, len = length(df))
    if (length(Pub) < length(df)) Pub <- rep(Pub, len = length(df))
    
    out <- numeric(max(length(zb), length(lb), length(ub), length(Plb), length(Pub)))
    out[ind] <- crps.norm(zb[ind], 0, 1, lb[ind], ub[ind], Plb[ind], Pub[ind])
    out[!ind] <- crps.t(zb[!ind], df[!ind], 0, 1, lb[!ind], ub[!ind], Plb[!ind], Pub[!ind])
    return(res + scale * out)
  }
  
  if (ind1 & ind2) {
    if (any(Plb + Pub > 1)){
      stop("Sum of 'lmass' and 'umass' exceeds 1.")
    }
    a <- ifelse(1 - Plb - Pub < 1e-12, 0,
                (1 - Plb - Pub) / (pt(ub, df) - pt(lb, df)))
    a[a > 1e12] <- NaN

    out_l <- sign(lb) * a^2 * sqrt(df) / (df - 1) *
      beta(0.5, df - 0.5) / beta(0.5, 0.5 * df)^2 *
      pbeta(df / (df + lb^2), df - 0.5, 0.5, lower.tail = FALSE) +
      ifelse(is.finite(lb),
             -lb * Plb^2 -
               2 * a * df / (df - 1) * (1 + lb^2/df) * dt(lb, df) * Plb,
             0)
    out_u <- -sign(ub) * a^2 * sqrt(df) / (df - 1) *
      beta(0.5, df - 0.5) / beta(0.5, 0.5 * df)^2 *
      pbeta(df / (df + ub^2), df - 0.5, 0.5, lower.tail = FALSE) +
      ifelse(is.finite(ub),
             ub * Pub^2 -
               2 * a * df / (df - 1) * (1 + ub^2/df) * dt(ub, df) * Pub,
             0)
    out_y <- zb * (2 * (a * (pt(zb, df) - pt(lb, df)) + Plb) - 1) +
      2 * a * df / (df - 1) * (1 + zb^2/df) * dt(zb, df)
  } else if (ind1 & !ind2) {
    a <- ifelse(1 - Plb < 1e-12, 0, (1 - Plb) / (1 - pt(lb, df)))
    a[a > 1e12] <- NaN
    
    out_l <- sign(lb) * a^2 * sqrt(df) / (df - 1) *
      beta(0.5, df - 0.5) / beta(0.5, 0.5 * df)^2 *
      pbeta(df / (df + lb^2), df - 0.5, 0.5, lower.tail = FALSE) +
      ifelse(is.finite(lb),
             -lb * Plb^2 -
               2 * a * df / (df - 1) * (1 + lb^2/df) * dt(lb, df) * Plb,
             0)
    out_u <- -a^2 * sqrt(df) / (df - 1) *
      beta(0.5, df - 0.5) / beta(0.5, 0.5 * df)^2
    out_y <- zb * (2 * (1 - a * pt(zb, df, lower.tail = FALSE)) - 1) +
      2 * a * df / (df - 1) * (1 + zb^2/df) * dt(zb, df)
  } else if (!ind1 & ind2) {
    a <- ifelse(1 - Pub < 1e-12, 0, (1 - Pub) / pt(ub, df))
    a[a > 1e12] <- NaN
    
    out_l <- -a^2 * sqrt(df) / (df - 1) *
      beta(0.5, df - 0.5) / beta(0.5, 0.5 * df)^2
    out_u <- -sign(ub) * a^2 * sqrt(df) / (df - 1) *
      beta(0.5, df - 0.5) / beta(0.5, 0.5 * df)^2 *
      pbeta(df / (df + ub^2), df - 0.5, 0.5, lower.tail = FALSE) +
      ifelse(is.finite(ub),
             ub * Pub^2 -
               2 * a * df / (df - 1) * (1 + ub^2/df) * dt(ub, df) * Pub,
           0)
    out_y <- zb * (2 * a * pt(zb, df) - 1) +
      2 * a * df / (df - 1) * (1 + zb^2/df) * dt(zb, df)
  }
  
  return(res + scale * (out_y + out_l + out_u))
}
