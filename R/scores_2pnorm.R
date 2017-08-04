#' Calculating scores for the two-piece-normal distribution
#'
#' @param y vector of observations.
#' @param scale1,scale2 vectors of positive scale parameters.
#' @param location vector of location parameters.
#' @return A vector of score values.
#' @name scores_2pnorm
NULL

#' @rdname scores_2pnorm
#' @export
crps_2pnorm <- function(y, scale1, scale2, location = 0) {
  if (!identical(location, 0)) y <- y - location
  z <- y
  z[z < 0] <- 0
  y[y > 0] <- 0
  s <- scale1 + scale2
  s[s == 0] <- Inf
  crps_gtcnorm(y, scale = scale1, upper = 0, umass = scale2 / s) +
    crps_gtcnorm(z, scale = scale2, lower = 0, lmass = scale1 / s)
}

#' @rdname scores_2pnorm
#' @export
logs_2pnorm <- function(y, scale1, scale2, location = 0) {
  -log(f2pnorm(y, location, scale1, scale2))
}
  


check_crps_2pnorm <- function(input) {
  required <- c("y", "location", "scale1", "scale2")
  checkNames1(required, names(input))
  checkNumeric(input)
  checkVector(input)
  
  if (any(input$scale1 <= 0))
    stop("Parameter 'scale1' contains non-positive values.")
  if (any(input$scale2 <= 0))
    stop("Parameter 'scale2' contains non-positive values.")
}

check_logs_2pnorm <- check_crps_2pnorm
