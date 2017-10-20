#' Calculating scores for the uniform distribution
#'
#' @param y vector of observations.
#' @param min,max lower and upper limits of the distribution. Must be finite.
#' @param lmass,umass vectors of point masses in \code{min} and \code{max}
#'  respectively.
#' @return A vector of score values.
#' @name scores_unif
NULL


#' @rdname scores_unif
#' @export
crps_unif <- function(y, min = 0, max = 1, lmass = 0, umass = 0) {
  if (identical(min, 0) && identical(max, 1)) {
    z <- y
    z[z < 0] <- 0
    z[z > 1] <- 1
    a <- 1 - (lmass + umass)
    a[a < 0] <- NaN
    abs(y - z) + z^2 * a - z * (1 - 2 * lmass) +
      a^2 / 3 + (1 - lmass) * umass
  } else {
    min[!is.finite(min)] <- NaN
    max[!is.finite(max)] <- NaN
    scale <- max - min
    scale[scale < 0] <- NaN
    if (all(scale > 0, na.rm = TRUE)) {
      scale * crps_unif((y - min) / scale, lmass = lmass, umass = umass)
    } else {
      out <- scale * crps_unif((y - min) / scale, lmass = lmass, umass = umass)
      ind <- scale == 0
      out[ind] <- rep_len(abs(y - min), length(out))[ind]
      out
    }
  }
}

#' @rdname scores_unif
#' @export
logs_unif <- function(y, min = 0, max = 1) {
  -dunif(y, min, max, log = TRUE)
}

#' @rdname scores_unif
#' @export
dss_unif <- function(y, min = 0, max = 1) {
  min[min > max] <- NaN
  m <- 0.5 * (min + max)
  s <- (max - min) / sqrt(12)
  ((y - m) / s)^2 + 2*log(s)
}
  


check_crps_unif <- function(input) {
  required <- c("y", "min", "max", "lmass", "umass")
  checkNames1(required, names(input))
  checkNumeric(input)
  checkVector(input)
  
  if (any(input$min > input$max))
    stop("Parameter 'min' contains greater values than parameter 'max'.")
  if (any(input$lmass < 0 | input$lmass > 1))
    stop("Parameter 'lmass' contains values not in [0, 1].")
  if (any(input$umass < 0 | input$umass > 1))
    stop("Parameter 'umass' contains values not in [0, 1].")
  if (any(input$lmass + input$umass > 1))
    stop("Values in 'lmass' and 'umass' add up to more than 1.")
}

check_logs_unif <- function(input) {
  required <- c("y", "min", "max")
  checkNames1(required, names(input))
  checkNumeric(input)
  checkVector(input)
  
  if (any(input$min > input$max))
    stop("Parameter 'min' contains greater values than parameter 'max'.")
}
