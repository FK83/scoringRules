#' Calculating scores for the log-logistic distribution
#'
#' @param y vector of observations.
#' @param locationlog vector of location parameters on the log scale.
#' @param scalelog vector of positive scale parameters on the log scale.
#' @return A vector of score values.
#' @name scores_llogis
NULL

#' @rdname scores_llogis
#' @export
crps_llogis <- function(y, locationlog, scalelog) {
  scalelog[scalelog < 0 | scalelog >= 1] <- NaN
  y1 <- y
  y1[y1 < 0] <- 0
  p <- plogis(log(y1), locationlog, scalelog)
  c1 <- y * (2 * p - 1)
  c2 <- 2 * exp(locationlog) * beta(1 + scalelog, 1 - scalelog)
  c3 <- (1 - scalelog) / 2 - pbeta(p, 1 + scalelog, 1 - scalelog)
  c1 + c2 * c3
}

#' @rdname scores_llogis
#' @export
logs_llogis <- function(y, locationlog, scalelog)
  -log(fllogis(y, locationlog, scalelog))
