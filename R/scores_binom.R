#' Calculating scores for the binomial distribution
#'
#' @param y vector of observations.
#' @inheritParams stats::pbinom
#' @return A vector of score values.
#' @name scores_binom
#' @importFrom stats pbinom dbinom
NULL

#' @rdname scores_binom
#' @export
crps_binom <- function(y, size, prob) {
  n_param <- max(length(size), length(prob))
  n_y <- length(y)
  size <- rep(size, length.out = n_param)
  prob <- rep(prob, length.out = n_param)
  
  if (n_y <= n_param) {
    y <- rep(y, length.out = n_param)
    sapply(
      seq_along(y),
      function(i) {
        y <- y[i]
        size <- size[i]
        prob <- prob[i]
        if (anyNA(c(y, size, prob))) return(y * size * prob)
        size_rounded <- round(size)
        tol <- .Machine$double.eps^0.5
        size <- if (abs(size - size_rounded) < tol) {
          size_rounded
        } else {
          warning(sprintf("non-integer n = %.6f", size))
          return(NaN)
        }
        if (size >= 0) {
          x <- seq.int(0, size, 1)
          w <- dbinom(x, size, prob)
          a <- pbinom(x, size, prob) - 0.5 * w
          2 * sum(w * ((y < x) - a) * (x - y))
        } else {
          NaN
        }
      }
    )
  } else {
    list_param <- lapply(
      seq_along(size),
      function(i) {
        size <- size[i]
        prob <- prob[i]
        if (anyNA(c(size, prob))) {
          typeNA <- size * prob
          return(list(x = typeNA, w = typeNA, a = typeNA))
        }
        size_rounded <- round(size)
        tol <- .Machine$double.eps^0.5
        listNaN <- list(x = NaN, w = NaN, a = NaN)
        size <- if (abs(size - size_rounded) < tol) {
          size_rounded
        } else {
          warning(sprintf("non-integer n = %.6f", size))
          return(listNaN)
        }
        if (size >= 0) {
          x <- seq.int(0, size, 1)
          w <- dbinom(x, size, prob)
          a <- pbinom(x, size, prob) - 0.5 * w
          list(x = x, w = w, a = a)
        } else {
          listNaN
        }
      }
    )
    list_param <- rep(list_param, length.out = n_y)
    sapply(
      seq_along(y),
      function(i) {
        with(list_param[[i]], 2 * sum(w * ((y[i] < x) - a) * (x - y[i])))
      }
    )
  }
}


#' @rdname scores_binom
#' @export
logs_binom <- function(y, size, prob) {
  -dbinom(y, size, prob, log = TRUE)
}

check_crps_binom <- function(input) {
  required <- c("y", "size", "prob")
  checkNames1(required, names(input))
  checkNumeric(input)
  checkVector(input)
  
  if (any(input$size <= 0))
    stop("Parameter 'size' contains non-positive values.")
  if (any(input$prob > 1 | input$prob < 0))
    stop("Parameter 'prob' contains values not in [0, 1].")
}

check_logs_binom <- check_crps_binom
