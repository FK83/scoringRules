#' Calculating the CRPS for the Laplace distribution
#'
#' @param y vector of observations.
#' @param location vector of location parameters.
#' @param scale vector of positive scale parameters.
#' @return A vector of CRPS values.
#' @export
crps_lapl <- function(y, location = 0, scale = 1) {
  if (!identical(location, 0)) y <- y - location
  if (identical(scale, 1)) {
    abs_y <- abs(y)
    abs_y + exp(-abs_y) - 0.75
  } else {
    scale[scale < 0] <- NaN
    if (all(scale > 0, na.rm = TRUE)) {
      scale * crps_lapl(y / scale)
    } else {
      out <- scale * crps_lapl(y / scale)
      ind <- scale == 0
      out[ind] <- rep_len(abs(y), length(out))[ind]
      out
    }
  }
}
