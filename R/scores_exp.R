#' Calculating scores for the exponential distribution
#'
#' Calculating scores (CRPS, logarithmic score) for the exponential distribution, and the exponential
#' distribution with location-scale transformation and point mass in
#' \code{location}.
#'
#' @param y vector of observations.
#' @param rate vector of rates.
#' @param location vector of location parameters.
#' @param scale vector of positive scale parameters.
#' @param mass vector of point masses in \code{location}.
#' @return A vector of score values.
#' @name scores_exp
#' @importFrom stats pexp dexp
NULL

#' @rdname scores_exp
#' @export
crps_exp <- function(y, rate = 1) {
  abs(y) - (2 * pexp(y, rate) - 0.5) / rate
}

#' @rdname scores_exp
#' @export
crps_expM <- function(y, location = 0, scale = 1, mass = 0) {
  if (!identical(location, 0)) y <- y - location
  mass[mass < 0 | mass > 1] <- NaN
  a <- 1 - mass
  abs(y) - scale * a * (2 * pexp(y, 1 / scale) - 0.5 * a)
}

#' @rdname scores_exp
#' @export
logs_exp <- function(y, rate = 1) {
  -dexp(y, rate, log = TRUE)
}

#' @rdname scores_exp
#' @export
logs_exp2 <- function(y, location = 0, scale = 1) {
  -dexp(y - location, 1 / scale, log = TRUE)
}


check_crps_exp <- function(input) {
  required <- c("y", "rate")
  checkNames1(required, names(input))
  checkNumeric(input)
  checkVector(input)
  
  if (any(input$rate <= 0))
    stop("Parameter 'rate' contains non-positive values.")
}

check_crps_expM <- function(input) {
  required <- c("y", "location", "scale", "mass")
  checkNames1(required, names(input))
  checkNumeric(input)
  checkVector(input)
  
  if (any(input$rate <= 0))
    stop("Parameter 'scale' contains non-positive values.")
  if (any(input$mass < 0 | input$mass > 1))
    stop("Parameter 'mass' contains values not in [0, 1].")
}

check_logs_exp <- check_crps_exp
  
check_logs_exp2 <- function(input) {
  required <- c("y", "location", "scale")
  checkNames1(required, names(input))
  checkNumeric(input)
  checkVector(input)
  
  if (any(input$scale <= 0))
    stop("Parameter 'scale' contains non-positive values.")
}
