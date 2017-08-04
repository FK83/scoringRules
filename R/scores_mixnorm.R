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
