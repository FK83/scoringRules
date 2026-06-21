#' Ranked Probability Score
#' @description Computes the Ranked Probability Score (RPS) for a given vector or matrix of probabilities.
#' @param y vector of realizations, taking integer values between 1 and K. For the RPS, outcomes have an ordinal interpretation only (see details).
#' @param x vector or matrix (depending on \code{y}; see details) of probabilities.
#' @param tol non-negative numeric value specifying the allowed deviation from 1 of the total (row-wise) probability in \code{x}.
#' @details 
#' The RPS interprets the outcome variable as ordinal. That is, different outcome values can be ranked (e.g., \code{y=1} is smaller than \code{y=2}), but their numerical difference has no meaningful interpretation. 
#' For simplicity, the outcome \code{y} is coded as an integer here, with \code{y = 1} indicating the smallest possible realization and \code{y = K} indicating the largest possible realization. 
#' If \code{y} is a vector of length n >= 2, \code{x} should be given as a matrix
#' with n rows and K columns. If \code{y} has length 1, then \code{x} may be a vector of length K.
#' @examples
#' # Example with three outcome categories (a single observation)
#' p <- c(.3, .2, .5)
#' y <- 2
#' rps_probs(y, p)
#' 
#' # Example with three outcome categories (n = 2 observations)
#' p <- matrix(c(.2, .4, .4, .3, .6, .1), nrow = 2, byrow = TRUE)
#' y <- c(2, 3)
#' rps_probs(y, p)
#' @references
#' Original proposal of the RPS 
#'
#' Epstein, E.S. (1969): `A Scoring System for Probability Forecasts of Ranked Categories', Journal of Applied Meteorology and Climatology 8, 985-987.
#'
#' Application example (see esp. Section 4 for comments on the RPS' ordinal interpretation)
#'
#' Krueger, F. and L. Pavlova (2024): `Quantifying Subjective
#' Uncertainty in Survey Expectations', International Journal of Forecasting 40, 796-810, \doi{10.1016/j.ijforecast.2023.06.001}.
#' @export
rps_probs <- function(y, x, tol = .Machine$double.eps){
  input <- list(y = y, x = x, tol = tol)
  checkNumeric(input, infinite_exception = "tol")
  
  if (is_scalar(y)) {
    check_sizes_VecMat(input, c(y = "scalar", x = "vector", tol = "scalar"))
    check_p0(y, x, tol)
    rps0(y = y, x = x)
  } else {
    check_sizes_VecMat(input, c(y = "vector", x = "matrix", tol = "scalar"))
    check_p2(input)
    sapply(seq_along(y),
           function(i) rps0(y = y[i], x = x[i, ]))
  }
}

rps0 <- function(y, x){
  Px <- cumsum(x) # non-exceedance probabilities by category
  Py <- as.numeric(y <= seq_along(x)) # non-exceedance events by category
  sum((Px - Py)^2) # sum of Brier scores over categories
}

check_p0 <- function(y, x, tol){
  msg <- character(0)
  K <- length(x)
  if (K < 2){
    msg <- c(msg, "Number of outcome categories K must be >= 2.")
  }
  if (isFALSE(y %in% 1:K)){
    msg <- c(msg, "Unexpected input for 'y' (should be integer between 1 and K, where K is the number of categories).")
  }
  if (isTRUE(any(x < 0))) {
    msg <- c(msg, "Probabilities in 'x' must be nonnegative.")
  }
  if (isTRUE(any(tol < 0))) {
    msg <- c(msg, "Parameter 'tol' must be nonnegative.")
  }
  if (isTRUE(abs(sum(x) - 1) > tol)){
    msg <- c(msg, "Probabilities in 'x' should sum to one (increase 'tol' if a deviation is intended).")
  }
  if (length(msg) > 0){
    stop(paste(msg, collapse = "\n"))
  }
}

check_p2 <- function(input){
  sapply(seq_along(input$y),
         function(i) check_p0(y = input$y[i], 
                              x = input$x[i, ],
                              tol = input$tol))
}
