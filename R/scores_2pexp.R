#' Calculating scores for the two-piece-exponential distribution
#'
#' @param y vector of observations.
#' @param scale1,scale2 vectors of positive scale parameters.
#' @param location vector of location parameters.
#' @return A vector of score values.
#' @name scores_2pexp
NULL

#' @rdname scores_2pexp
#' @export 
crps_2pexp <- function(y, scale1, scale2, location = 0) {
  if (!identical(location, 0)) y <- y - location
  z <- y
  z[z < 0] <- 0
  y[y > 0] <- 0
  s <- scale1 + scale2
  s[s == 0] <- Inf
  crps_expM(-y, scale = scale1, mass = scale2 / s) +
    crps_expM(z, scale = scale2, mass = scale1 / s)
}

#' @rdname scores_2pexp
#' @export 
logs_2pexp <- function(y, scale1, scale2, location = 0)
  -log(f2pexp(y, location, scale1, scale2))
