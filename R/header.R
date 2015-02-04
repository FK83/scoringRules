# !!! Important !!! Set orientation = 1 for positive orientation, to - 1 for negative orientation
# (only used in functions that appear in the package: xx.parametric and xx.sample, where xx = (ls, crps, dss, qs)
orientation <- -1

################################################################################
### xx.sample

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

ls.sample <- function(dat, y, bw = NULL){
  message("Using the log score with kernel density estimation tends to be fragile -- see KLTG (2015) for details.")
  if (is.null(bw)) bw <- bw.SJ(dat)
  out <- lsmixnC(w = rep(1/length(dat), length(dat)), m = dat, s = rep(bw, length(dat)), y = y)
  return(orientation*out)
}

dss.sample <- function(dat, y){
  message("Note: The DS Score considers only the first two moments of the simulated distribution.")
  out <- dss(mean(dat), sd(dat), y)
  return(orientation*out)
}

################################################################################
### parametric

# Header function (CRPS based on parametric family)
crps.parametric <- function(y, family, ...){
  ind <- c(match(family, names(list_inputChecks)),
           match(paste("check.", family, sep=""), list_inputChecks)
  )
  ind <- ind[!is.na(ind)]
  ind <- unique(ind)
  if (length(ind) > 1) {
    stop("Ambiguous choice of parametric family - see details section of ?crps.parametric for a list of available choices.")
  } else if (length(ind) == 0) {
    stop("Could not find parametric family - see details section of ?crps.parametric for a list of available choices.")  
  }
  
  input <- list(y = y, ...)
  inputCheck <- eval(parse(text = list_inputChecks[ind]))
  input <- inputCheck(input)
    
  ind <- names(list_inputChecks)[ind]
  crps <- eval(parse(text = list_crpsFunctions[ind]))
  formals(crps) <- input
  out <- crps()
  
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
