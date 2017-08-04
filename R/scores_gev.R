#' Calculating scores for the generalized extreme value distribution
#'
#' @param y vector of observations.
#' @param shape vector of positive shape parameters.
#' @param location vector of location parameters.
#' @param scale vector of positive scale parameters.
#' @return A vector of score values.
#' @name scores_gev
NULL

#' @rdname scores_gev
#' @export
crps_gev <- function(y, shape, location = 0, scale = 1) {
  shape[shape >= 1] <- NaN
  if (!identical(location, 0)) y <- y - location
  if (!identical(scale, 1)) {
    scale[scale < 0] <- NaN
    y <- y / scale
  } 
  
  if (any(ind <- abs(shape) < 1e-12, na.rm = TRUE)) {
    if (length(y) < length(shape)) y <- rep_len(y, length(shape))
    out <- rep_len(NaN, length(y))
    
    out[ind] <- -y[ind] - digamma(1) - log(2) - 2 *
      if (requireNamespace("gsl", quietly = TRUE)) {
        gsl::expint_Ei(-exp(-y[ind]))
      } else {
        warning(paste("The exponential integral is approximated using the 'integrate' function.",
                      "Consider installing the 'gsl' package to leverage a more accurate implementation.",
                      sep = "\n"))
        sapply(-exp(-y[ind]), function(upper) {
          integrate(function(x) exp(x)/x, -Inf, upper)$value
        })
      }
    out[!ind] <- crps_gev(y[!ind], shape[!ind])
  } else {
    x <- 1 + shape * y
    x[x < 0] <- 0
    x <- x^(-1/shape)
    c1 <- 2 * exp(-x) - 1
    out <- (y + 1/shape) * c1 + gamma(1 - shape) / shape *
      (2 * pgamma(x, 1 - shape) - 2^shape)
  }
  return(scale * out)
}

#' @rdname scores_gev
#' @export
logs_gev <- function(y, shape, location = 0, scale = 1) {
  -log(fgev(y, location, scale, shape))
}
  


check_crps_gev <- function(input) {
  required <- c("y", "location", "scale", "shape")
  checkNames1(required, names(input))
  checkNumeric(input)
  checkVector(input)
  
  if (any(input$scale <= 0))
    stop("Parameter 'scale' contains non-positive values.")
}

check_logs_gev <- check_crps_gev
