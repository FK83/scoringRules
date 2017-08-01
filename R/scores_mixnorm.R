#' Calculating the CRPS for a mixture of normal distribution.
#'
#' @param y vector of observations.
#' @param m matrix of mean parameters.
#' @param s matrix of scale parameters.
#' @param w matrix of weights.
#' @param exact if \code{TRUE} calculates the analytical solution; otherwise
#'  numerical integration is used.
#' @param rel_tol only used if \code{exact} is \code{FALSE}; relative tolerance
#'  for numerical integration.
#' @return A vector of CRPS values.
#' @export
crps_mixnorm = function(y, m, s, w, exact = TRUE, rel_tol = 1e-6){
  if (exact == TRUE){
    out <- sapply(seq_along(y), function(i) crpsmixnC(w[i, ], m[i, ], s[i, ], y[i]))
  } else {
    out <- sapply(seq_along(y), function(i) crps.mixnorm.int(y[i], m[i, ], s[i, ], w[i, ], rel_tol))
  }
  return(out)
}
