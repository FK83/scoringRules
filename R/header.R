# !!! Important !!! Set orientation = 1 for positive orientation, to - 1 for negative orientation
# (only used in functions that appear in the package: xx and xx_sample, where xx = (crps, logs, qs)
orientation <- -1

#' @export crps logs crps_sample logs_sample flapl fnorm f2pexp fmixnorm f2pnorm ft fllapl flogis fllogis fexp fgev fgpd
#' @importFrom Rcpp evalCpp
#' @importFrom methods existsFunction
#' @importFrom stats bw.nrd dbeta dexp dgamma dlnorm dlogis dnbinom dnorm dpois dt dunif integrate
#' @importFrom stats pbeta pexp pgamma plnorm plogis pnbinom pnorm ppois pt punif
#' @useDynLib scoringRules

################################################################################
### xx_sample

crps_sample <- function(y, dat, method = "edf", w = NULL, bw = NULL, 
                        num_int = FALSE, show_messages = TRUE) {
    input <- list(y = y, dat = dat)
    if (!is.null(w)) {
      input$w <- w
    }
    if (!is.null(bw)) {
      input$bw <- bw
    }
    check.sample(input)
    
    # Throw error if weights are not equal
    if (method == "kde" & !is.null(w)) {
      warning("Observation weights not implemented for kernel density estimation - edf used instead.")
      method <- "edf"
    }
    
    if (method == "edf") {
      out <- crps.edf(dat = dat, y = y, w = w)
    } else if (method == "kde") {
      # Bandwidth
      if (is.null(bw)) {
        bw <- bw.nrd(dat)
      }
      if (num_int == FALSE) {
        # Exact formula
        out <- crps.kdens(dat = dat, y = y, bw = bw)
      } else {
        # Numerical integration
        FF <- function(x) {
          sapply(x, function(zz) mean(pnorm(zz, mean = dat, sd = bw)))
        }
        aux1 <- integrate(function(s) FF(s)^2,-Inf, y, rel.tol = 1e-6)$value
        aux2 <- integrate(function(s) (1 - FF(s))^2, y, Inf, rel.tol = 1e-6)$value
        out <- aux1 + aux2
        # Message
        if (show_messages) 
          message("Used numerical integration for CRPS computation (tolerance = 1e-6).")
      }
    } else {
      stop("Unexpected choice of method - please select either edf or kde")
    }
    return(out)
  }

logs_sample <- function(y, dat, bw = NULL, show_messages = TRUE) {
  input <- list(y = y, dat = dat)
  if (!is.null(bw)) {
    input$bw <- bw
  }
  check.sample(input)
  
  if (show_messages)
    message("Using the log score with kernel density estimation tends to be fragile -- see KLTG (2016) for details.")
  if (is.null(bw)) bw <- bw.nrd(dat)
  w <- rep(1 / length(dat), length(dat))
  s <- rep(bw, length(dat))
  out <- lsmixnC(w = w, m = dat, s = s, y = y)
  return(out)
}

################################################################################
### parametric

crps <- function(y, family, ...) {
  family <- checkFamily(family, "crps")
  checkInput <- get(paste0("check.", family))
  calculateCRPS <- get(paste0("crps.", family))
  
  if (is.list(y)) {
    input <- c(y, list(...))
  } else {
    input <- list(y = y, ...)
  }
  
  input <- checkInput(input)
  out <- do.call(calculateCRPS, input)
  
  if (any(is.na(out))) {
    warning("Missing CRPS values. Probably due to numerical instabilities as a result of extreme parameter choices.")
  } else if (any(out < 0)) {
    warning("Negative CRPS values. Check parameter combinations and/or contact package maintainer(s).")
  }
  
  return(out)
}

logs <- function(y, family, ...) {
  family <- checkFamily(family, "ls")
  checkInput <- get(paste0("check.", family))
  calculateLS <- get(paste0("ls.", family))
  
  if (is.list(y)) {
    input <- c(y, list(...))
  } else {
    input <- list(y = y, ...)
  }
  
  input <- checkInput(input)
  out <- do.call(calculateLS, input)
  
  return(out)
}

################################################################################