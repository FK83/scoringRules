#' Calculating scores for the Poisson distribution
#'
#' @param y vector of observations.
#' @inheritParams stats::ppois
#' @return A vector of score values.
#' @name scores_pois
#' @importFrom stats ppois dpois
NULL

#' @rdname scores_pois
#' @export
crps_pois <- function(y, lambda) {
  c1 <- (y - lambda) * (2*ppois(y, lambda) - 1)
  c2 <- 2*dpois(floor(y), lambda) -
    exp(-2*lambda) * (besselI(2*lambda, 0) + besselI(2*lambda, 1))
  return(c1 + lambda*c2)
}

#' @rdname scores_pois
#' @export
logs_pois <- function(y, lambda) {
  -dpois(y, lambda, log = TRUE)
}

#' @rdname scores_pois
#' @export
dss_pois <- function(y, lambda) {
  lambda[lambda <= 0] <- NaN
  (y - lambda)^2 / lambda + log(lambda)
}


check_crps_pois <- function(input) {
  required <- c("y", "lambda")
  checkNames1(required, names(input))
  checkNumeric(input)
  checkVector(input)
  
  if (any(input$lambda <= 0))
    stop("Parameter 'lambda' contains non-positive values.")
}

check_logs_pois <- check_crps_pois
