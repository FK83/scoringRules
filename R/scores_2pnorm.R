#' Calculating the CRPS for the two-piece-normal distribution
#'
#' @param y vector of observations.
#' @param scale1,scale2 vectors of positive scale parameters.
#' @param location vector of location parameters.
#' @return A vector of CRPS values.
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
