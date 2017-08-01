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
logs_unif <- function(y, min, max)
  -dunif(y, min, max, log=TRUE)
