#' Calculating the CRPS for the generalized extreme value distribution
#'
#' @param y vector of observations.
#' @param shape vector of positive shape parameters.
#' @param location vector of location parameters.
#' @param scale vector of positive scale parameters.
#' @return A vector of CRPS values.
#' @export
crps_gev <- function(y, shape, location = 0, scale = 1) {
  shape[shape >= 1] <- NaN
  if (!identical(location, 0)) y <- y - location
  if (!identical(scale, 1)) y <- y / (scale[scale < 0] <- NaN)
  
  if (any(ind <- abs(shape) < 1e-12, na.rm = TRUE)) {
    if (length(z) < length(shape)) z <- rep_len(z, length(shape))
    out <- rep_len(NaN, length(z))
    
    out[ind] <- -z[ind] - digamma(1) - log(2) - 2 *
      if (requireNamespace("gsl", quietly = TRUE)) {
        gsl::expint_Ei(-exp(-z[ind]))
      } else {
        warning(paste("The exponential integral is approximated using the 'integrate' function.",
                      "Consider installing the 'gsl' package to leverage a more accurate implementation.",
                      sep = "\n"))
        sapply(-exp(-z[ind]), function(upper) {
          integrate(function(x) exp(x)/x, -Inf, upper)$value
        })
      }
    out[!ind] <- crps_gev(z[!ind], shape[!ind])
  } else {
    x <- 1 + shape * z
    x[x < 0] <- 0
    x <- x^(-1/shape)
    c1 <- 2 * exp(-x) - 1
    out <- (z + 1/shape) * c1 + gamma(1 - shape) / shape *
      (2 * pgamma(x, 1 - shape) - 2^shape)
  }
  return(scale * out)
}
