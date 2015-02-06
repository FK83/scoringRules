# !!! Important !!! Set orientation = 1 for positive orientation, to - 1 for negative orientation
# (only used in functions that appear in the package: xx.parametric and xx.sample, where xx = (ls, crps, dss, qs)
orientation <- -1

################################################################################
### xx.sample

crps.sample <- function(y, dat, method = "edf", w = NULL, bw = NULL){
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

qs.sample <- function(y, dat, bw = NULL){
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

ls.sample <- function(y, dat, bw = NULL){
  message("Using the log score with kernel density estimation tends to be fragile -- see KLTG (2015) for details.")
  if (is.null(bw)) bw <- bw.SJ(dat)
  out <- lsmixnC(w = rep(1/length(dat), length(dat)), m = dat, s = rep(bw, length(dat)), y = y)
  return(orientation*out)
}

dss.sample <- function(y, dat){
  message("Note: The DS Score considers only the first two moments of the simulated distribution.")
  out <- dnorm(y, mean(dat), sd(dat), log=TRUE)
  return(orientation*out)
}

################################################################################
### parametric

crps.parametric <- function(y, family, ...){
  ind <- c(match(family, names(list_crpsFunctions)),
           match(paste("crps.", family, sep=""), list_crpsFunctions)
  )
  ind <- ind[!is.na(ind)]
  ind <- unique(ind)
  ind <- names(list_crpsFunctions)[ind]
  if (length(ind) > 1) {
    stop("Ambiguous choice of parametric family - see details section of ?crps.parametric for a list of available choices.")
  } else if (length(ind) == 0) {
    stop("Could not find parametric family - see details section of ?crps.parametric for a list of available choices.")  
  }
  
  input <- list(y = y, ...)
  inputCheck <- eval(parse(text = list_inputChecks[ind]))
  input <- inputCheck(input)
  
  crps <- eval(parse(text = list_crpsFunctions[ind]))
  formals(crps) <- input
  out <- crps()
  
  return(orientation*out)
}

qs.parametric <- function(family, parameters, y){
  ind <- c(match(family, names(list_qsFunctions)),
           match(paste("qs.", family, sep=""), list_qsFunctions)
  )
  ind <- ind[!is.na(ind)]
  ind <- unique(ind)
  ind <- names(list_qsFunctions)[ind]
  if (length(ind) > 1) {
    stop("Ambiguous choice of parametric family - see details section of ?qs.parametric for a list of available choices.")
  } else if (length(ind) == 0) {
    stop("Could not find parametric family - see details section of ?qs.parametric for a list of available choices.")  
  }
  
  input <- list(y = y, ...)
  inputCheck <- eval(parse(text = list_inputChecks[ind]))
  input <- inputCheck(input)
  
  qs <- eval(parse(text = list_qsFunctions[ind]))
  formals(qs) <- input
  out <- qs()
  
  return(orientation*out)
}

ls.parametric <- function(y, family, ...){
  ind <- c(match(family, names(list_lsFunctions)),
           match(paste("ls.", family, sep=""), list_lsFunctions)
  )
  ind <- ind[!is.na(ind)]
  ind <- unique(ind)
  ind <- names(list_lsFunctions)[ind]
  if (length(ind) > 1) {
    stop("Ambiguous choice of parametric family - see details section of ?ls.parametric for a list of available choices.")
  } else if (length(ind) == 0) {
    stop("Could not find parametric family - see details section of ?ls.parametric for a list of available choices.")  
  }
  
  input <- list(y = y, ...)
  inputCheck <- eval(parse(text = list_inputChecks[ind]))
  input <- inputCheck(input)
  
  logscore <- eval(parse(text = list_lsFunctions[ind]))
  formals(logscore) <- input
  out <- logscore()
  
  return(orientation*out)
}

dss.parametric <- function(y, family, ...) {
  ind <- c(match(family, names(list_dssFunctions)),
           match(paste("dss.", family, sep=""), list_dssFunctions)
  )
  ind <- ind[!is.na(ind)]
  ind <- unique(ind)
  ind <- names(list_dssFunctions)[ind]
  if (length(ind) > 1) {
    stop("Ambiguous choice of parametric family - see details section of ?dss.parametric for a list of available choices.")
  } else if (length(ind) == 0) {
    stop("Could not find parametric family - see details section of ?dss.parametric for a list of available choices.")  
  }
  
  message("Note: The DS Score considers only the first two moments of the chosen distribution.")
  
  input <- list(y = y, ...)
  inputCheck <- eval(parse(text = list_inputChecks[ind]))
  input <- inputCheck(input)
  
  dss <- eval(parse(text = list_dssFunctions[ind]))
  formals(dss) <- input
  out <- dss()
  
  return(orientation*out)
}

################################################################################

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
