#' Calculating scores for the gamma distribution
#'
#' @param y vector of observations.
#' @param shape vector of positive shape parameters.
#' @param rate an alternative way to specify the scale.
#' @param scale vector of positive scale parameters.
#' @return A vector of score values.
#' @name scores_gamma
#' @importFrom stats pgamma dgamma
NULL

#' @rdname scores_gamma
#' @export
crps_gamma <- function(y, shape, rate = 1, scale = 1/rate) {
  # code fragments from stats::pgamma
  if (!missing(rate) && !missing(scale)) {
    
    if (all(abs(rate * scale - 1) < 1e-15) )
      warning("specify 'rate' or 'scale' but not both")
    else stop("specify 'rate' or 'scale' but not both")
  }
  p1 <- pgamma(y, shape, scale = scale)
  p2 <- pgamma(y, shape + 1, scale = scale)
  y * (2 * p1 - 1) - scale * (shape * (2 * p2 - 1) + 1 / beta(.5, shape))
}

#' @rdname scores_gamma
#' @export
logs_gamma <- function(y, shape, rate = 1, scale = 1/rate) {
  # code copied from stats::dgamma
  if (!missing(rate) && !missing(scale)) {
    if (abs(rate * scale - 1) < 1e-15) 
      warning("specify 'rate' or 'scale' but not both")
    else stop("specify 'rate' or 'scale' but not both")
  }
  -dgamma(y, shape, scale = scale, log = TRUE)
}
  
