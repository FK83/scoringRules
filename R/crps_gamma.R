#' Calculating the CRPS for the gamma distribution
#'
#' @param y vector of observations.
#' @param shape vector of positive shape parameters.
#' @param rate an alternative way to specify the scale.
#' @param scale vector of positive scale parameters.
#' @return A vector of CRPS values.
#' @export
crps_gamma <- function(y, shape, rate = 1, scale = 1/rate) {
  if (!missing(rate) && !missing(scale)) { # from stats::pgamma
    if (all(abs(rate * scale - 1) < 1e-15) )
      warning("specify 'rate' or 'scale' but not both")
    else stop("specify 'rate' or 'scale' but not both")
  }
  p1 <- pgamma(y, shape, scale = scale)
  p2 <- pgamma(y, shape + 1, scale = scale)
  y * (2 * p1 - 1) - scale * (shape * (2 * p2 - 1) + 1 / beta(.5, shape))
}
