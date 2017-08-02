################################################################################
#' @export crps logs crps_sample logs_sample flapl fnorm f2pexp fmixnorm f2pnorm ft fllapl flogis fllogis fexp fgev fgpd es_sample vs_sample
#' @importFrom Rcpp evalCpp
#' @importFrom methods existsFunction
#' @importFrom stats bw.nrd dbeta dexp dgamma dlnorm dlogis dnbinom dnorm dpois dt dunif integrate
#' @importFrom stats pbeta pexp pgamma plnorm plogis pnbinom pnorm ppois pt punif
#' @useDynLib scoringRules, .registration = TRUE
################################################################################
### parametric

crps <- function(y, ...) UseMethod("crps")

logs <- function(y, ...) UseMethod("logs")

#' @export
#' @export crps.numeric
crps.numeric <- function(y, family, ...) {
  family <- checkFamily(family, "crps")
  checkInput <- get(paste0("check.", family))
  calculateCRPS <- get(paste0("crps.", family))
  
  input <- list(y = y, ...)
  input <- checkInput(input)
  out <- do.call(calculateCRPS, input)
  
  if (any(is.na(out))) {
    warning("Missing CRPS values. Probably due to numerical instabilities as a result of extreme parameter choices.")
  } else if (any(out < 0)) {
    warning("Negative CRPS values. Check parameter combinations and/or contact package maintainer(s).")
  }
  
  return(out)
}

#' @export
#' @export logs.numeric
logs.numeric <- function(y, family, ...) {
  family <- checkFamily(family, "ls")
  checkInput <- get(paste0("check.", family))
  calculateLS <- get(paste0("ls.", family))
  
  input <- list(y = y, ...)
  input <- checkInput(input)
  out <- do.call(calculateLS, input)
  
  return(out)
}

################################################################################
