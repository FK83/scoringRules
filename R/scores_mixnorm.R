#' Calculating scores for a mixture of normal distributions.
#'
#' @param y vector of observations.
#' @param m matrix of mean parameters (rows represent observations, columns represent mixture components).
#' @param s matrix of scale parameters (same structure as for \code{m}).
#' @param w matrix of weights (same structure as for \code{m}; row sums must equal one).
#' @param exact if \code{TRUE} calculates the analytical solution; otherwise
#'  numerical integration is used.
#' @param rel_tol only used if \code{exact} is \code{FALSE}; relative tolerance
#'  for numerical integration.
#' @return A vector of score values.
#' @name scores_mixnorm
#' @examples
#' 
#' # Example 1: 100 observations, 15 mixture components
#' mval <- matrix(rnorm(100*15), nrow = 100)
#' sdval <- matrix(rgamma(100*15, shape = 2), nrow = 100)
#' weights <- matrix(rep(1/15, 100*15), nrow = 100)
#' y <- rnorm(100)
#' crps1 <- crps_mixnorm(y = y, m = mval, s = sdval, w = weights)
#' crps2 <- crps_mixnorm(y = y, m = mval, s = sdval, w = weights, exact = FALSE)
#' 
#' # Example 2: 2 observations, 10000 mixture components
#' mval <- matrix(rnorm(2*10000), nrow = 2)
#' sdval <- matrix(rgamma(2*10000, shape = 2), nrow = 2)
#' weights <- matrix(rep(1/10000, 2*10000), nrow = 2)
#' y <- rnorm(2)
#' # With many mixture components, non-exact evaluation is much faster
#' system.time(crps1 <- crps_mixnorm(y = y, m = mval, s = sdval, w = weights))
#' system.time(crps2 <- crps_mixnorm(y = y, m = mval, s = sdval, w = weights, exact = FALSE))
NULL

#' @rdname scores_mixnorm
#' @export
crps_mixnorm = function(y, m, s, w, exact = TRUE, rel_tol = 1e-6){
  if (exact == TRUE){
    sapply(seq_along(y),
           function(i) crpsmixnC(w[i, ], m[i, ], s[i, ], y[i]))
  } else {
    sapply(seq_along(y),
           function(i) crps_mixnorm_int(y[i], m[i, ], s[i, ], w[i, ], rel_tol))
  }
}

crps_mixnorm_int <- function(y, m, s, w, rel_tol){
  Fmix <- function(z){
    sapply(z, function(r) sum(w*pnorm((r-m)/s)))
  }
  crps_int(y, Fmix, -Inf, Inf, rel_tol)
}

#' @rdname scores_mixnorm
#' @export
logs_mixnorm <- function(y, m, s, w) {
  sapply(seq_along(y), function(i) lsmixnC(w[i, ], m[i, ], s[i, ], y[i]))
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
