#' Calculating scores for the negative binomial distribution
#'
#' @param y vector of observations.
#' @inheritParams stats::pnbinom
#' @return A vector of score values.
#' @name scores_nbinom
#' @importFrom stats pnbinom dnbinom
NULL

#' @rdname scores_nbinom
#' @export
crps_nbinom <- function(y, size, prob) {
  if (!requireNamespace("hypergeo", quietly = TRUE)) {
    stop(paste(
      "Calculations require an implementation of the gaussian hypergeometric function.",
      "Please install the following package: hypergeo (>= 1.0)",
      sep = "\n"))
  }
  c1 <- y * (2 * pnbinom(y, size, prob) - 1)
  c2 <- (1 - prob) / prob ^ 2
  c3 <- prob * (2 * pnbinom(y - 1, size + 1, prob) - 1) + Re(hypergeo::hypergeo(size + 1, 0.5, 2,-4 * c2))
  return(c1 - size * c2 * c3)
}

#' @rdname scores_nbinom
#' @export
logs_nbinom <- function(y, size, prob) 
  -dnbinom(y, size, prob, log = TRUE)
