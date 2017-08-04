#' Calculating scores for the Laplace distribution
#'
#' @param y vector of observations.
#' @param location vector of location parameters.
#' @param scale vector of positive scale parameters.
#' @return A vector of score values.
#' @name scores_lapl
NULL

#' @rdname scores_lapl
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

#' @rdname scores_lapl
#' @export
logs_lapl <- function(y, location = 0, scale = 1) {
  -log(flapl(y, location, scale))
}


check_crps_lapl <- function(input) {
  required <- c("y", "location", "scale")
  checkNames1(required, names(input))
  checkNumeric(input)
  checkVector(input)
  
  if (any(input$scale <= 0))
    stop("Parameter 'scale' contains non-positive values.")
}

check_logs_lapl <- check_crps_lapl
