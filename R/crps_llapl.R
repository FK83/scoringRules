#' Calculating the CRPS for the log-Laplace distribution
#'
#' @param y vector of observations.
#' @param locationlog vector of location parameters on the log scale.
#' @param scalelog vector of positive scale parameters on the log scale.
#' @return A vector of CRPS values.
#' @export
crps_llapl <- function(y, locationlog, scalelog) {
  scalelog[scalelog < 0 | scalelog >= 1] <- NaN
  y1 <- y
  y1[y1 < 0] <- 0
  z <- (log(y1) - locationlog) / scalelog
  p <- 0.5 + 0.5 * sign(z) * pexp(abs(z))
  c1 <- y*(2*p - 1)
  c2 <- ifelse (z < 0,
                (1 - (2*p)^(1 + scalelog)) / (1 + scalelog),
                - (1 - (2*(1-p))^(1 - scalelog)) / (1 - scalelog)
  )
  c3 <- scalelog / (4 - scalelog^2)
  c1 + exp(locationlog) * (c2 + c3)
}
