#' Calculating scores for the beta distribution
#'
#' @param y vector of observations.
#' @param shape1,shape2 vectors of positive shape parameters.
#' @param lower,upper vectors of lower and upper limits of the distribution. Must be finite.
#' @return A vector of score values.
#' @importFrom stats pbeta dbeta
#' @name scores_beta
NULL

#' @rdname scores_beta
#' @export
crps_beta <- function(y, shape1, shape2, lower = 0, upper = 1) {
  if (identical(lower, 0) && identical(upper, 1)) {
    z <- y
    z[z < 0] <- 0
    z[z > 1] <- 1
    c1 <- y * (2*pbeta(y, shape1, shape2) - 1)
    c2 <- shape1 / (shape1 + shape2)
    c3 <- 1 - 2 * pbeta(y, shape1 + 1, shape2)
    c4 <- 2 / shape1 * beta(2 * shape1, 2 * shape2) / beta(shape1, shape2)^2
    if (any(ind <- !is.finite(c4))) {
      c4[ind] <- sqrt(shape2 / (pi * shape1 * (shape1 + shape2)))[ind]  # stirling's approximation
    }
    c1 + c2 * (c3 - c4)
  } else {
    lower[!is.finite(lower)] <- NaN
    upper[!is.finite(upper)] <- NaN
    scale <- upper - lower
    scale[scale < 0] <- NaN
    if (all(scale > 0, na.rm = TRUE)) {
      scale * crps_beta((y - lower) / scale, shape1, shape2)
    } else {
      out <- scale * crps_beta((y - lower) / scale, shape1, shape2)
      ind <- scale == 0
      out[ind] <- rep_len(abs(y - lower), length(out))[ind]
      out
    }
  }
}

#' @rdname scores_beta
#' @export
logs_beta <- function(y, shape1, shape2, lower = 0, upper = 1) {
  if (identical(lower, 0) & identical(upper, 1)) {
    -dbeta(y, shape1, shape2, log = TRUE)
  } else {
    lower[!is.finite(lower)] <- NaN
    upper[!is.finite(upper)] <- NaN
    scale <- upper - lower
    scale[scale <= 0] <- NaN
    -dbeta((y - lower) / scale, shape1, shape2, log = TRUE) + log(scale)
  }
}

#' @rdname scores_beta
#' @export
dss_beta <- function(y, shape1, shape2, lower = 0, upper = 1) {
  if (!identical(lower, 0)) y <- y - lower
  shape1[shape1 <= 0] <- NaN
  shape2[shape2 <= 0] <- NaN
  scale <- upper - lower
  scale[scale <= 0] <- NaN
  r <- shape2 / shape1
  m <- scale / (1 + r)
  v <- m^2 * r / (1 + shape1 + shape2)
  (y - m)^2 / v + log(v)
}


check_crps_beta <- function(input) {
  required <- c("y", "shape1", "shape2", "lower", "upper")
  checkNames1(required, names(input))
  checkNumeric(input)
  checkVector(input)
  
  if (any(input$shape1 <= 0))
    stop("Parameter 'shape1' contains non-positive values.")
  if (any(input$shape1 <= 0))
    stop("Parameter 'shape1' contains non-positive values.")
  if (any(input$lower > input$upper))
    stop("Parameter 'lower' contains greater values than parameter 'upper'.")
}

check_logs_beta <- check_crps_beta
