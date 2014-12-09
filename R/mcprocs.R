# !!! Important !!! Set orientation = 1 for positive orientation, to - 1 for negative orientation
# (only used in functions that appear in the package: xx.parametric and xx.sample, where xx = (ls, crps, dss, qs)
orientation <- -1

# Simulate Data from the Fox/West Model
sim.fw <- function(a, s, n, ndraws){
  # Compute parameters of marginal distribution (inverse Gamma)
  beta.ig <- c(0.5*(n+2), 0.5*n*s)
  # Compute parameters of corresponding t distribution
  beta.t <- c(beta.ig[2]/beta.ig[1], 2*beta.ig[1]) # variance, degrees of freedom
  # Simulate data
  sig2 <- rep(0, ndraws)
  sig2[1] <- 1
  psi <- 1/rgamma(ndraws, 0.5*(n+3), 0.5*n*s*(1-a^2))
  v <- a + sqrt( psi/(n*s) ) * rnorm(ndraws)
  for (j in 2:ndraws){
    sig2[j] <- psi[j] + (v[j]^2) * sig2[j-1]
  }  
  # Return list
  return(list(par.ig = beta.ig, sim = sig2))
}

# Make the ergodic dist for our simulation design
ergdist <- function(s, n){
  # Make t distribution parameters
  s <- sqrt(n*s/(n + 2))
  df <- n + 2
  # Compute properties of ergodic dist
  pdf <- function(z){
    sapply(z, function(l) dt(l/s, df)/s)
  }  
  cdf <- function(z){
    sapply(z, function(l) pt(l/s, df))
  } 
  list(pdf = pdf, cdf = cdf, s = s, df = df)
}

# Numerical integration for CRPS
crps.int <- function(F, y, rel.tol = 1e-6){
  s1 <- integrate(function(s) F(s)^2, -Inf, y, rel.tol = rel.tol, subdivisions = 1000L, stop.on.error = TRUE)$value
  s2 <- integrate(function(s) (1-F(s))^2, y, Inf, rel.tol = rel.tol, subdivisions = 1000L, stop.on.error = TRUE)$value
  return(- (s1 + s2))
}

# CRPS of a (weighted) empirical distribution
crps.edf <- function(dat, y, w = NULL){
  n <- length(dat)
  # Set uniform weights unless specified otherwise
  if (is.null(w)){
    x <- .Internal(sort(dat, FALSE))
    out <- sapply(y, function(s) 2 / n^2 * sum((n * (s < x) - 1:n + 0.5) * (x - s)))
  } else {
    if (length(w) != n) stop()
    ord <- order(dat)
    x <- dat[ord]
    w <- w[ord]
    p <- c(0, cumsum(w))
    out <- sapply(y, function(s) 2 * sum((w * (s < x) - 0.5 * (p[2:(n+1)]^2 - p[1:n]^2)) * (x - s)))
  }
  return(-out)
}

# CRPS for a mixture of normals
crps.mixn = function (m, s, y, w = NULL, exact = TRUE, rel.tol = 1e-6){
  n <- length(m)
  if (is.null(w)) 
    w <- rep(1/n, n)
  if (exact == TRUE){
    return(crpsmixnC(w, m, s, y))
  } else {
    Fmix = function(z){
      sapply(z, function(r) sum(w*pnorm((r-m)/s)))
    }
    return(crps.int(Fmix, y, rel.tol = 1e-6))
  }
}

# CRPS for a t distribution
crps.t <- function(m, s, df, y){
  z <- (y-m)/s
  c1 <- (2*s*dt(z, df)) * ((df + z^2)/(df - 1)) + z*s*(2*pt(z, df) - 1)
  c2 <- ( s*(4*sqrt(df))/(df - 1) ) * (beta(0.5, df - 0.5)/(beta(0.5, 0.5*df)^2))
  return(- c1 + 0.5*c2)
}

# CRPS for a GEV distribution
crps.gev <- function(loc, sc, sh, y){
  if (abs(sh) < 1e-12) {message("Small shape parameter. Numerical integration was used for CPRS (tolerance = 1e-6)")
                        pgev.tmp <- function(x){return(exp(-exp(-(x - loc)/sc)))}
                        return(-integrate(function(x) (pgev.tmp(x) - (y <= x))^2, lower=-Inf, upper=Inf, rel.tol = 1e-6, subdivisions = 1000L, stop.on.error = TRUE)$value)
  }
  else if (y < loc-sc/sh) {
    return(-((loc - sc/sh - y) - sc/sh*( (2^sh-2)*gamma(1-sh))))
  }
  else {
    logFy <- (1 + sh*(y - loc)/sc)^{-1/sh}
    return(-((loc - sc/sh - y)*(1 - 2*exp(-logFy)) - 2*sc/sh * gamma(1-sh)*(2^(sh-1) - pgamma(logFy,shape=1-sh))))}
}

# CRPS for a GPD distribution
crps.gpd <- function(loc, sc, sh, y){
  z <- (y-loc)/sc
  if (abs(sh) < 1e-12 & y > loc) {
    Fy <- Fgpd(y,loc=loc,sc=sc,sh=sh) 
    return(-(y-loc-sc*(2*Fy - 0.5)))
  }
  else if (abs(sh) < 1e-12 & y < loc) {(-1)*return((-1)*integrate(function(x) (Vectorize(Fgpd)(x,loc=loc,sc=sc,sh=sh) - (y <= x))^2, y, Inf, rel.tol=1e-6)$value )}
  else if (sh > 1e-12 & y < loc) {return((-1)*integrate(function(x) (Vectorize(Fgpd)(x,loc=loc,sc=sc,sh=sh) - (y <= x))^2, -Inf, Inf, rel.tol=1e-6)$value)}
  else if (sh == 1) {return((-1)*integrate(function(x) (Vectorize(Fgpd)(x,loc=loc,sc=sc,sh=sh) - (y <= x))^2, -Inf, Inf, rel.tol=1e-6)$value)}
  else if (sh < -1e-12 & y >= loc-sc/sh) {return((-1)*integrate(function(x) (Vectorize(Fgpd)(x,loc=loc,sc=sc,sh=sh) - (y <= x))^2, -Inf, y, rel.tol=1e-6)$value)}
  else if (sh < -1e-12 & y < loc) {return((-1)*integrate(function(x) (Vectorize(Fgpd)(x,loc=loc,sc=sc,sh=sh) - (y <= x))^2, y, Inf, rel.tol=1e-6)$value )}
  else {
    Fy <- Fgpd(y,loc=loc,sc=sc,sh=sh)
    return(-((y - loc + sc/sh)*(2*Fy - 1) - 2*sc/(sh*(sh-1))*(1/(sh-2) + (1-Fy)*(1+sh*z))))
  }
}

# CRPS for a truncated normal distribution
crps.tn <- function(m, s, lb, y) {
  z <- (y-m)/s; Phi.z <- pnorm(z)
  x.lb <- (lb-m)/s; Phi.lb <- 1-pnorm(x.lb)
  return(-(s/Phi.lb^2*(z*Phi.lb*(2*Phi.z+Phi.lb-2) + 2*dnorm(z)*Phi.lb - (1-pnorm(sqrt(2)*x.lb))/sqrt(pi))))
}

# CRPS for a log-normal distribution
crps.ln <- function(loc, sh, y){
  if (y < 0) {
    return((-1)*integrate(function(x) (plnorm(x,meanlog=loc,sdlog=sh) - (y <= x))^2, y, Inf, rel.tol=1e-6, subdivisions=1000L)$value)
  }
  else {
    z <- (log(y) - loc)/sh
    m <- exp(loc + sh^2/2)
    return(-(y*(2*pnorm(z)-1) - 2*m*(pnorm(z-sh) + pnorm(sh/sqrt(2))  - 1)))      
  }   
}

# CRPS via kernel density estimation
crps.kdens = function(dat, y, bw = NULL, exact = TRUE, rel.tol = 1e-6){
  n <- length(dat)
  if (is.null(bw)) {
    s <- rep(bw.SJ(dat), n)
  }
  else {
    s <- rep(bw, n)
  }
  return(crps.mixn(dat, s, y, exact = exact, rel.tol = 1e-6))
}

# Log Score for Mixture of Normals
ls.mixn <- function(m, s, y, w = NULL) {
  n <- length(m)
  if (is.null(w)){
    w <- rep(1/n, n)
  }
  return(lsmixnC(w, m, s, y))
}

# Log Score for Kernel Density Estimate
ls.kdens <- function(dat, y, bw = NULL) {
  if (is.null(bw)) 
    bw <- bw.SJ(dat)
  return(ls.mixn(dat, rep(bw, length(dat)), y))
}

##############################################
# qs.2pnorm - Compute QS of 2-piece normal predictive density
# Inputs:
# - rlz, scalar, realizing value
# - m,s1,s2, distribution parameters (see p2pnorm above)
# Output:
# - QS (the greater the better)
##############################################
qs.2pnorm <- function(m, s1, s2, y){
  2*d2pnorm(y, m, s1, s2) - (s1 + s2)/(sqrt(pi)*(s1 + s2)^2)
}

##############################################
# crps.2pnorm - Compute CRPS of 2-piece normal predictive density
# Same as qs.2pnorm above, but for CRPS
##############################################
crps.2pnorm <- function(m, s1, s2, y){
  aux1 <- function(y,m,s1,s2){
    (4*(s1^2)/(s1+s2)) * (  ((y-m)/s1)*pnorm((y-m)/s1)+dnorm((y-m)/s1)      )    
  }
  aux2 <- function(m,s1,s2){
    (2/(sqrt(pi)*(s1+s2)^2)) * ( sqrt(2)*s2*(s2^2-s1^2)-(s1^3+s2^3) )    
  }
  sc <- ifelse(y <= m,aux1(y,m,s1,s2)-(y-m)+aux2(m,s1,s2),aux1(y,m,s2,s1)+(y-m)*( ((s1-s2)^2-4*s2^2) / ((s1+s2)^2) )+aux2(m,s2,s1))
  return(-sc)
}

##############################################
# CDF and density functions of extreme value distributions, required for computation of LS, QS (GEV + GPD) and CRPS (GPD)
##############################################
# Density of a GEV distribution
fgev <- function(x, loc, sc, sh){
  z <- (x - loc)/sc 
  if (!(abs(sh)<1e-12) & sh*z <= -1) {
    return(0)
  }
  else if (abs(sh) < 1e-12) {
    tmp <- log(1/sc) - z - exp(-z)
    return(exp(tmp))
  }
  else {
    tmp <- (1 + sh*z)
    return(exp(log(1/sc) - tmp^(-1/sh) - (1/sh + 1)*log(tmp)))
  }
}

# CDF of a GPD distribution
Fgpd <- function(x, loc, sc, sh){
  z <- (x-loc)/sc 
  if (abs(sh) < 1e-12) {
    if(x < loc) {return(0)}
    else {return(1-exp(-z))} 
  }
  else {
    if (sh > 0 & x > loc) {return(1 - (1 + sh*z)^(-1/sh))}
    else if (sh < 0 & x > loc & x < loc - sc/sh) {return(1 - (1 + sh*z)^(-1/sh))}
    else if (sh < 0 & x > loc & x > loc - sc/sh) {return(1)}
    else {return(0)}
  }
}

# Density of a GPD distribution
fgpd <- function(x, loc, sc, sh){
  z <- (x-loc)/sc
  if (abs(sh) < 1e-12) {
    if(x < loc) {return(0)}
    else {return(1/sc*exp(-z))}
  }
  else {
    if (sh > 0 & x > loc) {return( 1/sc*(1+sh*z)^(-(sh+1)/sh) )}
    else if (sh < 0 & x > loc & x < loc - sc/sh) {return( 1/sc*(1+sh*z)^(-(sh+1)/sh))}
    else {return(0)}
  }
}

# Header function (CRPS based on sample of data)
crps.sample <- function(dat, y, method = "edf", w = NULL, bw = NULL){
	# Throw error if weights are not equal
	if (method == "kde" & !is.null(w)){
		warning("Observation weights not implemented for kernel density estimation - edf used instead.")
		method <- "edf"
	}

	if (method == "edf"){
		out <- crps.edf(dat = dat, y = y, w = w)
	} else if (method == "kde") {
		# Exact method for "small" samples
		xx <- (length(dat) < 10000)
		out <- crps.kdens(dat = dat, y = y, bw = bw, exact = xx)
		# Messages
		if (xx == FALSE) message("Large sample - used numerical integration for CPRS (tolerance = 1e-6).")
	} else {
		stop("Unexpected choice of method - please select either edf or kde")
	}
	
	return(orientation*out)
	
}

##############################################
# d2pnorm - compute density of 2-piece normal
# Input/output analogous to p2pnorm above, but for pdf
# Setting lg = TRUE yields log density
##############################################
d2pnorm <- function(x, m, s1, s2, lg=FALSE){
  a <- 2/(sqrt(2*pi)*(s1+s2))
  ret <- ifelse(x <= m,a*sqrt(2*pi)*s1*dnorm(x,m,s1),a*sqrt(2*pi)*s2*dnorm(x,m,s2))  
  if (lg==FALSE) ret else log(ret)
}

# Header function (CRPS based on parametric family)
crps.parametric <- function(family, parameters, y){
  input.check.parametric(family, parameters, y)
  if (family == "normal"){
    out <- crps.mixn(m = parameters$m, s = parameters$s, y = y)
  } else if (family == "mixture-normal"){
    xx <- (length(parameters$m) < 10000)
    if (xx == FALSE) message("Large sample - used numerical integration for crps computation (tolerance = 1e-6).")
    out <- crps.mixn(m = parameters$m, s = parameters$s, y = y, w = parameters$w, exact = xx)
  } else if (family == "t"){
    out <- crps.t(m = parameters$m, s = parameters$s, df = parameters$df, y)
  } else if (family == "two-piece-normal"){
    out <- crps.2pnorm(m = parameters$m, s1 = parameters$s1, s2 = parameters$s2, y = y)
  } else if (family == "truncated-normal"){
    out <- crps.tn(m = parameters$m, s = parameters$s, lb = parameters$lb, y = y)
  } else if (family == "log-normal"){
    out <- crps.ln(loc = parameters$loc, sh = parameters$sh, y = y)
  } else if (family == "gev"){
    out <- crps.gev(loc = parameters$loc, sc = parameters$sc, sh = parameters$sh, y = y)
  } else if (family == "gpd"){
    out <- crps.gpd(loc = parameters$loc, sc = parameters$sc, sh = parameters$sh, y = y)
  } 
  return(orientation*out)
}

qs.parametric <- function(family, parameters, y){
  input.check.parametric(family, parameters, y)
  if (family == "normal"){
    out <- qsmixnC(w = 1, m = parameters$m, s = parameters$s, y = y)
  } else if (family == "mixture-normal"){
    out <- qsmixnC(w = rep(1/length(parameters$m), length(parameters$m)), m = parameters$m, s = parameters$s, y = y)
  } else if (family == "t"){
    ff <- function(x){
      z <- (x-parameters$m)/parameters$s
      dt(z, parameters$df)/parameters$s
    }
    out <- 2*ff(y) - integrate(function(x) ff(x)^2, -Inf, Inf, rel.tol = 1e-6)$value 
  } else if (family == "two-piece-normal"){
    out <- qs.2pnorm(m = parameters$m, s1 = parameters$s1, s2 = parameters$s2, y = y)
  } else if (family == "truncated-normal"){
    ff <- function(x){
      if (x < parameters$lb) {return(0)}
      else {return(dnorm(x,parameters$m,parameters$s)/(1-pnorm(parameters$lb,parameters$m,parameters$s)))} 
    }
    out <- 2*ff(y) - integrate(Vectorize(function(x) ff(x)^2), -Inf, Inf, rel.tol = 1e-6)$value    
  } else if (family == "log-normal"){
    out <- 2*dlnorm(y,meanlog=parameters$loc,sdlog=parameters$sh) - integrate(function(x) dlnorm(x,meanlog=parameters$loc,sdlog=parameters$sh)^2, -Inf, Inf, subdivisions = 1000L, rel.tol = 1e-6)$value    
  } else if (family == "gev"){
    out <- 2*fgev(y, loc=parameters$loc, sc=parameters$sc, sh=parameters$sh) - integrate(Vectorize(function(x) fgev(x, loc=parameters$loc, sc=parameters$sc, sh=parameters$sh)^2), -Inf, Inf, subdivisions = 1000L, rel.tol = 1e-6)$value    
  } else if (family == "gpd"){
    out <- 2*fgpd(y, loc=parameters$loc, sc=parameters$sc, sh=parameters$sh) - integrate(Vectorize(function(x) fgpd(x, loc=parameters$loc, sc=parameters$sc, sh=parameters$sh)^2), -Inf, Inf, subdivisions = 1000L, rel.tol = 1e-6)$value    
  } 
  return(orientation*out)
}

qs.sample <- function(dat, y, bw = NULL){
	if (is.null(bw)) bw <- bw.SJ(dat)
	if (length(dat) < 10000){
		out <- qsmixnC(w = rep(1/length(dat), length(dat)), m = dat, s = rep(bw, length(dat)), y = y)
	} else {
		ff <- function(x){
			sapply(x, function(zz) mean(dnorm(zz, mean = dat, sd = bw)))
		}
		message("Large sample - used numerical integration for qs computation (tolerance = 1e-6).")
		out <- 2*ff(y) - integrate(function(x) ff(x)^2, -Inf, Inf, rel.tol = 1e-6)$value
	}
	return(orientation*out)
}

ls.parametric <- function(family, parameters, y){
  input.check.parametric(family, parameters, y)
  if (family == "normal"){
    out <- lsmixnC(w = 1, m = parameters$m, s = parameters$s, y = y)
  } else if (family == "mixture-normal"){
    out <- lsmixnC(w = rep(1/length(parameters$m), length(parameters$m)), m = parameters$m, s = parameters$s, y = y)
  } else if (family == "t"){
    ff <- function(x){
      z <- (x-parameters$m)/parameters$s
      dt(z, parameters$df)/parameters$s
    }
    out <- log(ff(y))
  } else if (family == "two-piece-normal"){
    out <- d2pnorm(x = y, m = parameters$m, s1 = parameters$s1, s2 = parameters$s2, lg = TRUE)
  } else if (family == "truncated-normal"){
    ff <- function(x){
      if (x < parameters$lb) {return(0)}
      else {return(dnorm(x,parameters$m,parameters$s)/(1-pnorm(parameters$lb,parameters$m,parameters$s)))} 
    }
    out <- log(ff(y))
  } else if (family == "log-normal"){
    out <- log(dlnorm(y,meanlog=parameters$loc,sdlog=parameters$sh))
  } else if (family == "gev"){
    out <- log(fgev(y, loc=parameters$loc, sc=parameters$sc, sh=parameters$sh)) 
  } else if (family == "gpd"){
    out <- log(fgpd(y, loc=parameters$loc, sc=parameters$sc, sh=parameters$sh))  
  } 
  return(orientation*out)
}

ls.sample <- function(dat, y, bw = NULL){
	message("Using the log score with kernel density estimation tends to be fragile -- see KLTG (2015) for details.")
	if (is.null(bw)) bw <- bw.SJ(dat)
	out <- lsmixnC(w = rep(1/length(dat), length(dat)), m = dat, s = rep(bw, length(dat)), y = y)
	return(orientation*out)
}

dss <- function(m, s, y) return(dnorm(y, mean = m, sd = s, log = TRUE))

dss.sample <- function(dat, y){
	message("Note: The DS Score considers only the first two moments of the simulated distribution.")
	out <- dss(mean(dat), sd(dat), y)
	return(orientation*out)
}

dss.parametric <- function(family, parameters, y){
  input.check.parametric(family, parameters, y)
  if (family == "normal"){
    out <- dss(m = parameters$m, s = parameters$s, y = y)
  } else if (family == "mixture-normal"){
    message("Note: The DS Score considers only the first two moments of the chosen distribution.")
    m <- sum(parameters$m * parameters$w)
    v <- sum( ( (parameters$s^2) + (parameters$m - m)^2 ) * parameters$w )
    out <- dss(m = m, s = sqrt(v), y = y)
  } else if (family == "t"){
    message("Note: The DS Score considers only the first two moments of the chosen distribution.")
    s <- sqrt( (parameters$s^2) * (parameters$df/(parameters$df - 2)) )
    out <- dss(m = parameters$m, s = s, y = y)
  } else if (family == "two-piece-normal"){
    message("Note: The DS Score considers only the first two moments of the chosen distribution.")
    # Mean and variance, see eq. (A.2) and (A.3) in Box A of Wallis (2004)
    m <- parameters$m + sqrt(2/pi)*(parameters$s2 - parameters$s1)
    v <- (1-2/pi) * (parameters$s2 - parameters$s1)^2 + parameters$s1 * parameters$s2
    out <- dss(m = m, s = sqrt(v), y = y)
  } else if (family == "truncated-normal"){
    message("Note: The DS Score considers only the first two moments of the chosen distribution.")
    z <- (parameters$lb-parameters$m)/parameters$s
    m <- parameters$m + parameters$s*dnorm(z)/(1-pnorm(z))
    v <- parameters$s^2*(1 - dnorm(z)/(1-pnorm(z))*(dnorm(z)/(1-pnorm(z))-z))
    out <- dss(m = m, s = sqrt(v), y = y)
  } else if (family == "log-normal"){
    message("Note: The DS Score considers only the first two moments of the chosen distribution.")
    m <- exp(parameters$loc + parameters$sh^2/2)
    v <- exp(2*parameters$loc + parameters$sh^2)*(exp(parameters$sh^2)-1)
    out <- dss(m = m, s = sqrt(v), y = y)
  } else if (family == "gev"){
    message("Note: The DS Score considers only the first two moments of the chosen distribution.")
    if (parameters$sh < 1 & abs(parameters$sh) > 1e-12) {m <- parameters$loc + parameters$sc*(gamma(1-parameters$sh)-1)/parameters$sh} 
    else if (abs(parameters$sh) <= 1e-12) {m <- parameters$loc + parameters$sc*0.5772157}
    else if (parameters$sh >= 1) {m <- Inf; message("First moment does not exist for this shape parameter.")}
    if (parameters$sh < 0.5 & abs(parameters$sh) > 1e-12) {v <- parameters$sc^2*(gamma(1-2*parameters$sh) - gamma(1-parameters$sh)^2)/parameters$sh^2}
    else if (abs(parameters$sh) <= 1e-12) {v <- parameters$sc^2*(pi^2)/6}
    else if (parameters$sh >= 0.5) {v <- Inf; message("Second moment does not exist for this shape parameter.")}
    out <- dss(m = m, s = sqrt(v), y = y)
  } else if (family == "gpd"){
    message("Note: The DS Score considers only the first two moments of the chosen distribution.")
    if (parameters$sh < 1) {m <- parameters$loc + parameters$sc/(1-parameters$sh)} 
    else if (parameters$sh >= 1) {m <- Inf; message("First moment does not exist for this shape parameter.")}
    if (parameters$sh < 0.5) {v <- parameters$sc^2/((1-parameters$sh)^2*(1-2*parameters$sh))}
    else if (parameters$sh >= 0.5) {v <- Inf; message("Second moment does not exist for this shape parameter.")}
    out <- dss(m = m, s = sqrt(v), y = y)
  } 
  return(orientation*out)
}

# Add error checks for parameters (positivity etc)
input.check.parametric <- function(family, parameters, y){
  # Check 1: Length of parameter vectors
  if (family != "mixture-normal" & any(sapply(parameters, function(x) length(x)) != 1)) stop("All elements of list 'parameters' must be numeric vectors of length one - please see ?crps.parametric for details.")
  if (family == "mixture-normal" & (sd(sapply(parameters, function(x) length(x))) > 0) ) stop("Parameter vectors for mixture-normal must have the same length - please see ?crps.parametric for details.")
  
  # Check 2: Consistency of family and parameter list
  if (family == "normal"){
    if (any(sort(names(parameters)) != c("m", "s"))) stop("Unexpected input format -- Family 'normal' requires parameters 'm' and 's'. Please see ?crps.parametric for details.")
  } else if (family == "two-piece-normal"){
    if (any(sort(names(parameters)) != c("m", "s1", "s2"))) stop("Unexpected input format -- Family 'two-piece-normal' requires parameters 'm', 's1'  and 's2'. Please see ?crps.parametric for details.")
  } else if (family == "mixture-normal"){
    if (any(sort(names(parameters)) != c("m", "s", "w"))) stop("Unexpected input format -- Family 'mixture-normal' requires parameters 'm', 's' and 'w'. Please see ?crps.parametric for details.")
  } else if (family == "t"){
    if (any(sort(names(parameters)) != c("df", "m", "s"))) stop("Unexpected input format -- Family 't' requires parameters 'm', 's' and 'df'. Please see ?crps.parametric for details.")
  } else if (family == "gev"){
    if (any(sort(names(parameters)) != c("loc", "sc", "sh"))) stop("Unexpected input format -- Family 'gev' requires parameters 'loc', 'sc' and 'sh'. Please see ?crps.parametric for details.")
    if (parameters$sc <= 0) stop("Invalid scale parameter")
  } else if (family == "truncated-normal"){
    if (any(sort(names(parameters)) != c("lb", "m", "s"))) stop("Unexpected input format -- Family 'truncated-normal' requires parameters 'lb', 'm' and 's'. Please see ?crps.parametric for details.")
    if (parameters$s <= 0) stop("Invalid scale parameter")
  } else if (family == "log-normal"){
    if (any(sort(names(parameters)) != c("loc", "sh"))) stop("Unexpected input format -- Family 'log-normal' requires parameters 'loc' and 'sh'. Please see ?crps.parametric for details.")
    if (parameters$sh <= 0) stop("Invalid shape parameter")
  } else if (family == "gpd"){
    if (any(sort(names(parameters)) != c("loc", "sc", "sh"))) stop("Unexpected input format -- Family 'gpd' requires parameters 'loc', 'sc' and 'sh'. Please see ?crps.parametric for details.")
    if (parameters$sc <= 0) stop("Invalid scale parameter")
  } else {
    stop("Unexpected choice of parametric family - see details section of ?crps.parametric for a list of available choices")
  }
  
  # Check 3: Length of y
  if (length(y) != 1) stop("y must be of length one - a vectorized version is _not_ implemented to avoid ambiguities.")
}
