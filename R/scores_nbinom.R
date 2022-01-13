#' Calculating scores for the negative binomial distribution
#'
#' @param y vector of observations.
#' @inheritParams stats::pnbinom
#' @return A vector of score values.
#' @details The mean of the negative binomial distribution is given by \code{mu} = \code{size}*(1-\code{prob})/\code{prob}.
#' @name scores_nbinom
#' @importFrom stats pnbinom dnbinom
NULL

#' @rdname scores_nbinom
#' @export
crps_nbinom <- function(y, size, prob, mu) {
  # check from stats::pnbinom
  if (!missing(mu)) {
    if (!missing(prob))
      stop("specify 'prob' or 'mu' but not both")
    prob <- size / (size + mu)
  }
  if (!requireNamespace("hypergeo", quietly = TRUE)) {
    stop(paste(
      "Calculations require an implementation of the gaussian hypergeometric function.",
      "Please install the following package: hypergeo (>= 1.0)",
      sep = "\n"))
  }
  c1 <- y * (2 * pnbinom(y, size, prob) - 1)
  c2 <- (1 - prob) / prob ^ 2
  c3 <- (prob * (2 * pnbinom(y - 1, size + 1, prob) - 1)
         #+ Re(hypergeo::hypergeo(size + 1, 0.5, 2, -4 * c2)))
         + hypergeo_0.5_2(size + 1, -4 * c2))
  return(c1 - size * c2 * c3)
}

hypergeo_0.5_2 <- function(a, z) {
  if (!identical(length(a), 1L) || !identical(length(z), 1L)) {
    l <- max(length(a), length(z))
    a <- rep_len(a, l)
    z <- rep_len(z, l)
    return(sapply(seq_along(a), function(i) {hypergeo_0.5_2(a[i], z[i])}))
  }
  # 2F1(a, 0.5, 2, z)
  res <- NA
  if (a > 2.5 && abs(a - round(a)) < 1e-12) {
    res <- (1 - z)^(1.5 - a) *
             hypergeo::hypergeo_taylor(2 - a, 1.5, 2, z, maxiter = a - 2)
  }
  if (!is.finite(res) && z > -1) {
    res <- hypergeo::f15.3.4(a, 0.5, 2, z)
  }
  if (!is.finite(res)) {
    res <- Re(hypergeo::f15.3.8(a, 0.5, 2, z))
  }
  res
}

#' @rdname scores_nbinom
#' @export
logs_nbinom <- function(y, size, prob, mu) {
  if (!missing(prob)) {
    if (!missing(mu))
      stop("specify 'prob' or 'mu' but not both")
    -dnbinom(y, size, prob, log = TRUE)
  } else {
    -dnbinom(y, size, mu = mu, log = TRUE)
  }
}

#' @rdname scores_nbinom
#' @export
dss_nbinom <- function(y, size, prob, mu) {
  size[size <= 0] <- NaN
  if (!missing(prob)) {
    if (!missing(mu))
      stop("specify 'prob' or 'mu' but not both")
    prob[prob <= 0] <- NaN
    prob[prob >= 1] <- NaN
    mu <- size * (1 - prob) / prob
    v <- mu / prob
  } else {
    mu[mu < 0] <- NaN
    v <- mu * (1 + mu / size)
  }
  (y - mu)^2 / v + log(v)
}


check_crps_nbinom <- function(input) {
  required <- list(c("y", "size", "prob"),
                   c("y", "size", "mu"))
  checkNames2(required, names(input))
  checkNumeric(input)
  checkVector(input)
  
  if (any(input$size <= 0))
    stop("Parameter 'size' contains non-positive values.")
  if ("prob" %in% names(input)) {
    if (any(input$prob > 1 | input$prob <= 0))
      stop("Parameter 'prob' contains values not in (0, 1].")
  }
  if ("mu" %in% names(input)) {
    if (any(input$mu < 0))
      stop("Parameter 'mu' contains negative values.")
  }
}

check_logs_nbinom <- check_crps_nbinom
