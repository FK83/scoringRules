#' Calculating the CRPS for the Poisson distribution
#'
#' @param y vector of observations.
#' @inheritParams stats::ppois
#' @return A vector of CRPS values.
#' @export
# poisson
crps_pois <- function(y, lambda) {
  c1 <- (y - lambda) * (2*ppois(y, lambda) - 1)
  c2 <- 2*dpois(floor(y), lambda) -
    exp(-2*lambda) * (besselI(2*lambda, 0) + besselI(2*lambda, 1))
  return(c1 + lambda*c2)
}
