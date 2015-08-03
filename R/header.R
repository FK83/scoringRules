# !!! Important !!! Set orientation = 1 for positive orientation, to - 1 for negative orientation
# (only used in functions that appear in the package: xx.parametric and xx.sample, where xx = (ls, crps, dss, qs)
orientation <- -1

################################################################################
### xx.sample

crps_sample <- function(y, dat, method = "edf", w = NULL, bw = NULL, num_int = FALSE){
	# Throw error if weights are not equal
	if (method == "kde" & !is.null(w)){
		warning("Observation weights not implemented for kernel density estimation - edf used instead.")
		method <- "edf"
	}

	if (method == "edf"){
		out <- crps.edf(dat = dat, y = y, w = w)
	} else if (method == "kde") {
		# Bandwidth
		if (is.null(bw)) bw <- bw.nrd(dat)
		if (num_int == FALSE){
		  # Exact formula
		  out <- crps.kdens(dat = dat, y = y, bw = bw)
		} else {
		  # Numerical integration
		  FF <- function(x){
			sapply(x, function(zz) mean(pnorm(zz, mean = dat, sd = bw)))
		  }
		  aux1 <- integrate(function(s) FF(s)^2, -Inf, y, rel.tol = 1e-6)$value
		  aux2 <- integrate(function(s) (1-FF(s))^2, y, Inf, rel.tol = 1e-6)$value
		  out <- -(aux1 + aux2)
		  # Message
          message("Used numerical integration for CPRS (tolerance = 1e-6).")
		}
	} else {
		stop("Unexpected choice of method - please select either edf or kde")
	}
	return(orientation*out)
}

qs_sample <- function(y, dat, bw = NULL, num_int = FALSE){
  if (is.null(bw)) bw <- bw.nrd(dat)
  if (num_int == FALSE){
    out <- qsmixnC(w = rep(1/length(dat), length(dat)), m = dat, s = rep(bw, length(dat)), y = y)
  } else {
    ff <- function(x){
      sapply(x, function(zz) mean(dnorm(zz, mean = dat, sd = bw)))
    }
    message("Used numerical integration for qs computation (tolerance = 1e-6).")
    out <- 2*ff(y) - integrate(function(x) ff(x)^2, -Inf, Inf, rel.tol = 1e-6)$value
  }
  return(orientation*out)
}

logs_sample <- function(y, dat, bw = NULL){
  message("Using the log score with kernel density estimation tends to be fragile -- see KLTG (2015) for details.")
  if (is.null(bw)) bw <- bw.nrd(dat)
  out <- lsmixnC(w = rep(1/length(dat), length(dat)), m = dat, s = rep(bw, length(dat)), y = y)
  return(orientation*out)
}

################################################################################
### parametric

crps <- function(y, family, ...){
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

qs <- function(y, family, ...){
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

logs <- function(y, family, ...){
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

################################################################################