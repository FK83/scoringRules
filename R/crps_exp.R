#' Calculating the CRPS for the exponential distribution
#'
#' Calculating the CRPS for the exponential distribution, and the exponential
#' distribution with location-scale transformation and point mass in
#' \code{location}.
#'
#' @param y vector of observations.
#' @param rate vector of rates.
#' @param location vector of location parameters.
#' @param scale vector of positive scale parameters.
#' @param mass vector of point masses in \code{location}.
#' @return A vector of CRPS values.
#' @export
crps_exp <- function(y, rate = 1) {
  abs(y) - (2 * pexp(y, rate) - 0.5) / rate
}

#' @rdname crps_exp
#' @export
crps_expM <- function(y, location = 0, scale = 1, mass = 0) {
  if (!identical(location, 0)) y <- y - location
  mass[mass < 0 | mass > 1] <- NaN
  a <- 1 - mass
  abs(y) - scale * a * (2 * pexp(y, 1 / scale) - 0.5 * a)
}
