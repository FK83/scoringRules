#' Calculating the CRPS for the log-normal distribution
#'
#' @param y vector of observations.
#' @param meanlog an alternative way to specify \code{locationlog}.
#' @param sdlog an alternative way to specify \code{scalelog}.
#' @param locationlog vector of location parameters on the log scale.
#' @param scalelog vector of positive scale parameters on the log scale.
#' @return A vector of CRPS values.
#' @export
crps_lnorm <- function(y, meanlog = 0, sdlog = 1,
                       locationlog = meanlog, scalelog = sdlog) {
  if (!missing(meanlog) && !missing(locationlog)) {
    if (all(abs(meanlog - locationlog) < 1e-15))
      warning("specify 'meanlog' or 'locationlog' but not both")
    else stop("specify 'meanlog' or 'locationlog' but not both")
  }
  if (!missing(sdlog) && !missing(scalelog)) {
    if (all(abs(sdlog - scalelog) < 1e-15))
      warning("specify 'sdlog' or 'scalelog' but not both")
    else stop("specify 'sdlog' or 'scalelog' but not both")
  }
  c1 <- y * (2 * plnorm(y, meanlog, sdlog) - 1)
  c2 <- 2 * exp(meanlog + 0.5 * sdlog^2)
  c3 <- plnorm(y, meanlog + sdlog^2, sdlog) + pnorm(sdlog / sqrt(2)) - 1
  c1 - c2 * c3
}
