#' Calculating scores for a mixture of normal distributions.
#'
#' @param y vector of observations.
#' @param m matrix of mean parameters; rows represent observations, columns represent mixture components.
#' @param s matrix of scale parameters; same structure as \code{m}.
#' @param w optional; matrix of non-negative weights; same structure as \code{m}.
#' @param rel_tol relative accuracy for numerical integration.
#' @return A vector of score values.
#' @details \code{logs_mixnorm} and \code{crps_mixnorm} calculate scores via analytical formulas. \code{crps_mixnorm_int} uses numerical integration for the CRPS; this can be faster if there are many mixture components (i.e., if \code{m}, \code{s} and \code{w} have many columns). See examples below.
#' @name scores_mixnorm
#' @examples
#' 
#' # Example 1: 100 observations, 15 mixture components
#' mval <- matrix(rnorm(100*15), nrow = 100)
#' sdval <- matrix(rgamma(100*15, shape = 2), nrow = 100)
#' weights <- matrix(rep(1/15, 100*15), nrow = 100)
#' y <- rnorm(100)
#' crps1 <- crps_mixnorm(y = y, m = mval, s = sdval, w = weights)
#' crps2 <- crps_mixnorm_int(y = y, m = mval, s = sdval, w = weights)
#' 
#' \dontrun{
#' # Example 2: 2 observations, 10000 mixture components
#' mval <- matrix(rnorm(2*10000), nrow = 2)
#' sdval <- matrix(rgamma(2*10000, shape = 2), nrow = 2)
#' weights <- matrix(rep(1/10000, 2*10000), nrow = 2)
#' y <- rnorm(2)
#' # With many mixture components, numerical integration is much faster
#' system.time(crps1 <- crps_mixnorm(y = y, m = mval, s = sdval, w = weights))
#' system.time(crps2 <- crps_mixnorm_int(y = y, m = mval, s = sdval, w = weights))
#' }
NULL

#' @rdname scores_mixnorm
#' @export
crps_mixnorm <- function(y, m, s, w = NULL) {
  if (is.null(w)) {
    w <- rep(1, ncol(m))
    sapply(seq_along(y), function(i) crpsmixnC(w, m[i, ], s[i, ], y[i]))
  } else {
    sapply(seq_along(y), function(i) crpsmixnC(w[i, ], m[i, ], s[i, ], y[i]))
  }
}

#' @rdname scores_mixnorm
#' @export
crps_mixnorm_int <- function(y, m, s, w = NULL, rel_tol = 1e-6) {
  if (is.null(w)) {
    constructF <- function(i) {
      function(x) {
        sapply(x, function(z) mean(pnorm(z, m[i, ], s[i, ])))
      }
    }
  } else {
    w[w < 0] <- NaN
    constructF <- function(i) {
      W <- sum(w[i, ])
      function(x) {
        sapply(x, function(z) sum(w[i, ] * pnorm(z, m[i, ], s[i, ])) / W)
      }
    }
  }
  sapply(
    seq_along(y),
    function(i) crps_int(y[i], constructF(i), -Inf, Inf, rel_tol)
  )
}

#' @rdname scores_mixnorm
#' @export
logs_mixnorm <- function(y, m, s, w = NULL) {
  if (is.null(w)) {
    w <- rep(1, ncol(m))
    sapply(seq_along(y), function(i) lsmixnC(w, m[i, ], s[i, ], y[i]))
  } else {
    sapply(seq_along(y), function(i) lsmixnC(w[i, ], m[i, ], s[i, ], y[i]))
  }
}

#' @rdname scores_mixnorm
#' @export
dss_mixnorm <- function(y, m, s, w = NULL) {
  if (is.null(w)) {
    w <- rep(1, ncol(m))
    sapply(seq_along(y), function(i) dssmixnC(w, m[i, ], s[i, ], y[i]))
  } else {
    sapply(seq_along(y), function(i) dssmixnC(w[i, ], m[i, ], s[i, ], y[i]))
  }
}
  

check_crps_mixnorm <- function(input) {
  required <- c("y", "m", "s", "w")
  checkNames1(required, names(input))
  checkNumeric(input)
  checkMatrix(input)
  
  if (any(input$s <= 0))
    stop("Parameter 's' contains non-positive values.")
  if (any(input$w < 0 | input$w > 1))
    stop("Parameter 'w' contains values not in [0, 1].")
  if (!isTRUE(all.equal(apply(input$w, 1, sum), rep(1, dim(input$w)[1]))))
    stop("Parameter 'w' contains weighting schemes which do not sum up to 1.")
}

check_logs_mixnorm <- check_crps_mixnorm
