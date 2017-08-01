#' Calculating scores for a mixture of normal distributions.
#'
#' @param y vector of observations.
#' @param m matrix of mean parameters.
#' @param s matrix of scale parameters.
#' @param w matrix of weights.
#' @param exact if \code{TRUE} calculates the analytical solution; otherwise
#'  numerical integration is used.
#' @param rel_tol only used if \code{exact} is \code{FALSE}; relative tolerance
#'  for numerical integration.
#' @return A vector of score values.
#' @name scores_mixnorm
NULL

#' @rdname scores_mixnorm
#' @export
crps_mixnorm = function(y, m, s, w, exact = TRUE, rel_tol = 1e-6){
  if (exact == TRUE){
    sapply(seq_along(y),
           function(i) crpsmixnC(w[i, ], m[i, ], s[i, ], y[i]))
  } else {
    sapply(seq_along(y),
           function(i) crps.mixnorm.int(y[i], m[i, ], s[i, ], w[i, ], rel_tol))
  }
}

#' @rdname scores_mixnorm
#' @export
logs_mixnorm <- function(y, m, s, w) {
  sapply(seq_along(y), function(i) lsmixnC(w[i, ], m[i, ], s[i, ], y[i]))
}
  
