#' Calculating scores for the hypergeometric distribution
#'
#' @param y vector of observations / numbers of white balls drawn without replacement from an urn which contains both black and white balls.
#' @inheritParams stats::phyper
#' @return A vector of score values.
#' @name scores_hyper
#' @importFrom stats phyper dhyper
NULL

#' @rdname scores_hyper
#' @export
crps_hyper <- function(y, m, n, k) {
  n_param <- max(length(m), length(n), length(k))
  n_y <- length(y)
  m <- rep(m, length.out = n_param)
  n <- rep(n, length.out = n_param)
  k <- rep(k, length.out = n_param)
  
  if (n_y <= n_param) {
    y <- rep(y, length.out = n_param)
    sapply(
      seq_along(y),
      function(i) {
        y <- y[i]
        m <- m[i]
        n <- n[i]
        k <- k[i]
        if (anyNA(c(y, m, n, k))) return(y * m * n * k)
        m_rounded <- round(m)
        n_rounded <- round(n)
        k_rounded <- round(k)
        tol <- .Machine$double.eps^0.5
        m <- if (abs(m - m_rounded) < tol) m_rounded else return(NaN)
        n <- if (abs(n - n_rounded) < tol) n_rounded else return(NaN)
        k <- if (abs(k - k_rounded) < tol) k_rounded else return(NaN)
        if (m >= 0 && n >= 0 && k >= 0) {
          if (k > m + n) return(NaN)
          x <- seq.int(max(0, k - n), min(k, m), 1)
          w <- dhyper(x, m, n, k)
          a <- phyper(x, m, n, k) - 0.5 * w
          2 * sum(w * ((y < x) - a) * (x - y))
        } else {
          NaN
        }
      }
    )
  } else {
    list_param <- lapply(
      seq_along(m),
      function(i) {
        m <- m[i]
        n <- n[i]
        k <- k[i]
        if (anyNA(c(m, n, k))) {
          typeNA <- m * n * k
          return(list(x = typeNA, w = typeNA, a = typeNA))
        }
        m_rounded <- round(m)
        n_rounded <- round(n)
        k_rounded <- round(k)
        tol <- .Machine$double.eps^0.5
        listNaN <- list(x = NaN, w = NaN, a = NaN)
        m <- if (abs(m - m_rounded) < tol) m_rounded else return(listNaN)
        n <- if (abs(n - n_rounded) < tol) n_rounded else return(listNaN)
        k <- if (abs(k - k_rounded) < tol) k_rounded else return(listNaN)
        if (m >= 0 && n >= 0 && k >= 0) {
          x <- seq.int(max(0, k - n), min(k, m), 1)
          w <- dhyper(x, m, n, k)
          a <- phyper(x, m, n, k) - 0.5 * w
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


#' @rdname scores_hyper
#' @export
logs_hyper <- function(y, m, n, k) {
  -dhyper(y, m, n, k, log = TRUE)
}


check_crps_hyper <- function(input) {
  required <- c("y", "m", "n", "k")
  checkNames1(required, names(input))
  checkNumeric(input)
  checkVector(input)
  
  if (any(input$m <= 0))
    stop("Parameter 'm' contains non-positive values.")
  if (any(input$n <= 0))
    stop("Parameter 'n' contains non-positive values.")
  if (any(input$k <= 0))
    stop("Parameter 'k' contains non-positive values.")
  if (any(input$k > input$m + input$n))
    stop("Parameter 'k' contains values larger than the respective values in 'm + n'.")
  if (any(abs(input$m - round(input$m)) > sqrt(.Machine$double.eps)))
    stop("Parameter 'm' is not a whole number.")
  if (any(abs(input$n - round(input$n)) > sqrt(.Machine$double.eps)))
    stop("Parameter 'n' is not a whole number.")
  if (any(abs(input$k - round(input$k)) > sqrt(.Machine$double.eps)))
    stop("Parameter 'k' is not a whole number.")
}

check_logs_hyper <- check_crps_hyper
