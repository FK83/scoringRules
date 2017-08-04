#' Calculating scores for the generalized Pareto distribution
#'
#' @param y vector of observations.
#' @param shape vector of positive shape parameters.
#' @param location vector of location parameters.
#' @param scale vector of positive scale parameters.
#' @param mass vector of point masses in \code{location}.
#' @return A vector of score values.
#' @name scores_gpd
NULL

#' @rdname scores_gpd
#' @export
# generalized pareto distribution
crps_gpd <- function(y, shape, location = 0, scale = 1, mass = 0) {
  shape[shape >= 1] <- NaN
  if (!identical(location, 0)) y <- y - location
  if (!identical(scale, 1)) {
    scale[scale < 0] <- NaN
    z <- y / scale
  } else {
    z <- y
  }
  mass[mass < 0 | mass > 1] <- NaN
  
  x <- 1 + shape * z
  x[x < 0] <- 0
  x <- x^(-1/shape)
  if (any(ind <- abs(shape) < 1e-12, na.rm = TRUE)) {
    x <- ifelse(ind & seq_along(x), exp(-z), x)
  }
  x[x > 1] <- 1
  #p <- 1 - x
  a <- 1 - mass
  b <- 1 - shape
  
  abs(y) - scale * a * (2 / b * (1 - x^b) - a / (2 - shape))
}

#' @rdname scores_gpd
#' @export
logs_gpd <- function(y, shape, location = 0, scale = 1) {
  -fgpd(y, location, scale, shape, 0, log = TRUE)
}


check_crps_gpd <- function(input) {
  required <- c("y", "location", "scale", "shape", "mass")
  checkNames1(required, names(input))
  checkNumeric(input)
  checkVector(input)
  
  if (any(input$scale <= 0))
    stop("Parameter 'scale' contains non-positive values.")
  if (any(input$mass < 0 | input$mass > 1))
    stop("Parameter 'mass' contains values not in [0, 1].")
}

check_logs_gpd <- function(input) {
  required <- c("y", "location", "scale", "shape")
  checkNames1(required, names(input))
  checkNumeric(input)
  checkVector(input)
  
  if (any(input$scale <= 0))
    stop("Parameter 'scale' contains non-positive values.")
}
