#' Ranked Probability Score
#' @description Computes the Ranked Probability Score (RPS) for a given vector or matrix of probabilities.
#' @param y vector of realizations, taking integer values between 1 and K. For the RPS, outcomes have an ordinal interpretation only (see details).
#' @param x vector or matrix (depending on \code{y}; see details) of probabilities
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
#' Krüger, F., and L. Pavlova (2023): `Quantifying Subjective
#' Uncertainty in Survey Expectations', International Journal of Forecasting, forthcoming, \doi{j.ijforecast.2023.06.001}.
#' @export
rps_probs <- function(y, x){
  input <- list(y = y, dat = x)
  if (identical(length(y), 1L) && is.vector(x)) {
    check_p1(input)
    rps0(y = y, x = x)
  } else {
    check_p2(input)
    sapply(seq_along(y),
           function(i) rps0(y = y[i], x = x[i, ]))
  }
}

rps0 <- function(y, x){
  P <- cumsum(x)
  K <- length(x)
  if (y == 1){
    sum((1-P)^2)
  } else {
    ( sum(P[1:(y-1)]^2) + sum((1-P[y:K])^2) )
  }  
}

check_p0 <- function(y, x){
  msg1 <- msg2 <- msg3 <- ""
  K <- length(x)
  err <- 0
  if (K < 2){
    msg1 <- "Number of outcome categories K must be >= 2"
    err <- err + 1
  }
  if (!(y %in% 1:K)){
    msg2 <- "Unexpected input for y (should be integer between 1 and K, where K is the number of categories)"
    err <- err + 1
  } 
  if (any(x < 0) || (sum(x) > 1)){
    msg3 <- "Probabilities in x must be positive and sum to one"
    err <- err + 1
  }
  if (err > 0){
    stop(paste(msg1, msg2, msg3))
  }
}

check_p1 <- function(input){
  check_sample(input)
  check_p0(y = input$y, x = input$dat)
}

check_p2 <- function(input){
  check_sample2(input)
  sapply(seq_along(input$y),
         function(i) check_p0(y = input$y[i], 
                              x = input$dat[i, ]))
}