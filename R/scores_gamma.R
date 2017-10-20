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
  if (!missing(rate) && !missing(scale))
    stop("specify 'rate' or 'scale' but not both")
  p1 <- pgamma(y, shape, scale = scale)
  p2 <- pgamma(y, shape + 1, scale = scale)
  y * (2 * p1 - 1) - scale * (shape * (2 * p2 - 1) + 1 / beta(.5, shape))
}

#' @rdname scores_gamma
#' @export
logs_gamma <- function(y, shape, rate = 1, scale = 1/rate) {
  if (!missing(rate) && !missing(scale))
    stop("specify 'rate' or 'scale' but not both")
  -dgamma(y, shape, scale = scale, log = TRUE)
}

#' @rdname scores_gamma
#' @export
dss_gamma <- function(y, shape, rate = 1, scale = 1/rate) {
  if (!missing(rate) && !missing(scale))
    stop("specify 'rate' or 'scale' but not both")
  ms <- sqrt(shape)
  scale[scale <= 0] <- NaN
  s <- ms * scale
  (y / s - ms)^2 + 2*log(s)
}
# mean = shape * scale
# sd = sqrt(shape) * scale


check_crps_gamma <- function(input) {
  required <- list(c("y", "shape", "rate"),
                   c("y", "shape", "scale"))
  checkNames2(required, names(input))
  checkNumeric(input)
  checkVector(input)
  
  if (any(input$shape <= 0))
    stop("Parameter 'shape' contains non-positive values.")
  if ("rate" %in% names(input)) {
    if (any(input$rate <= 0))
      stop("Parameter 'rate' contains non-positive values.")
  }
  if ("scale" %in% names(input)) {
    if (any(input$scale <= 0))
      stop("Parameter 'scale' contains non-positive values.")
  }
}

check_logs_gamma <- check_crps_gamma
