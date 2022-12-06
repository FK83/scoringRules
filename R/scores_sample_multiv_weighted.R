#' Weighted Multivariate Scoring Rules for Simulated Forecast Distributions
#' 
#' Compute weighted versions of multivariate scores \eqn{S(y, dat)}, where \eqn{S} is a
#' proper scoring rule, \eqn{y} is a d-dimensional realization vector and 
#' \eqn{dat} is a simulated sample of multivariate forecasts. Three scores
#' are available: The energy score, a score based on a Gaussian kernel 
#' (\link{mmds_sample}, see details below) and the variogram score of order \eqn{p}.
#' 
#' @param y realized values (numeric vector of length d).
#' @param dat numeric matrix of data
#' (columns are simulation draws from multivariate forecast distribution).
#' @param a numeric vector of of length d containing lower bounds for the indicator 
#' weight function \code{w(z) = 1{a_{1} < z_{1} < b_{1}, ..., a_{d} < z_{d} < b_{d}}}.
#' @param b numeric vector of of length d containing lower bounds for the indicator 
#' weight function \code{w(z) = 1{a_{1} < z_{1} < b_{1}, ..., a_{d} < z_{d} < b_{d}}}.
#' @param w numeric vector of weights for forecast draws (length equal to number of columns of \code{dat})
#' @param w_vs numeric matrix of weights for \code{dat} used in the variogram
#' score. This matrix must be square and symmetric, with all elements being non-negative.
#' If no weights are specified, constant weights (with all elements of \code{w_vs} 
#' equal to one) are used.
#' @param p order of variogram score. Standard choices include \eqn{p = 1} and
#' \eqn{p = 0.5}.
#' 
#' 
#' @details
#' In the input matrix \code{dat} each column is expected to represent a sample
#' from the multivariate forecast distribution, the number of rows of \code{dat}
#' thus has to match the length of the observation vector \code{y}, and the
#' number of columns of \code{dat} is the number of simulated samples.
#' 
#' 
#' @return
#' Value of the score. \emph{A lower score indicates a better forecast.}
#' 
#' @references
#' \emph{Threshold-weighted scores}
#' 
#' Allen, S., Ginsbourger, D. and J. Ziegel (2022):
#' `Evaluating forecasts for high-impact events using transformed kernel scores',
#' arXiv preprint arXiv:2202.12732. 
#' \doi{10.48550/arXiv.2202.12732}
#' 
#' \emph{Outcome-weighted scores:}
#'  
#' Holzmann, H. and B. Klar (2017):
#' `Focusing on regions of interest in forecast evaluation',
#' \emph{Annals of Applied Statistics} 11, 2404-2431. 
#' \doi{10.1214/17-AOAS1088}
#' 
#' @author Sam Allen
#' 
#' 
#' @name scores_sample_multiv_weighted
NULL

################################################################################
# energy score
#' @rdname scores_sample_multiv_weighted
#' @export
twes_sample <- function(y, dat, a = -Inf, b = Inf, w = NULL) {
  if (length(a) == 1) {
    a <- rep(a, length(y))
  }else if (length(b) == 1) {
    b <- rep(b, length(y))
  }else if (length(a) != length(y) | length(b) != length(y)) {
    stop("The vector of weight values is not the same length as the vector of realizations")
  }
  v_y <- pmin(pmax(y, a), b)
  v_dat <- sapply(1:nrow(dat), function(m) pmin(pmax(dat[m, ], a), b))
  v_dat <- t(v_dat)
  score <- es_sample(y = v_y, dat = v_dat, w)
  return (score)
}

#' @rdname scores_sample_multiv_weighted
#' @export
owes_sample <- function(y, dat, a = -Inf, b = Inf, w = NULL) {
  weight_func <- function(x) as.numeric(x > a & x < b)
  w_y <- sapply(y, weight_func)
  w_dat <- array(sapply(dat, weight_func), dim(dat))
  score <- es_sample(y, dat, w = w_dat)*w_y 
  return (score)
}

################################################################################
# MMD score
#' @rdname scores_sample_multiv_weighted
#' @export
twmmds_sample <- function(y, dat, a = -Inf, b = Inf, w = NULL) {
  if (length(a) == 1) {
    a <- rep(a, length(y))
  }else if (length(b) == 1) {
    b <- rep(b, length(y))
  }else if (length(a) != length(y) | length(b) != length(y)) {
    stop("The vector of weight values is not the same length as the vector of realizations")
  }
  v_y <- pmin(pmax(y, a), b)
  v_dat <- sapply(1:nrow(dat), function(m) pmin(pmax(dat[m, ], a), b))
  v_dat <- t(v_dat)
  score <- mmds_sample(y = v_y, dat = v_dat, w)
  return (score)
}

#' @rdname scores_sample_multiv_weighted
#' @export
owmmds_sample <- function(y, dat, a = -Inf, b = Inf, w = NULL) {
  weight_func <- function(x) as.numeric(x > a & x < b)
  w_y <- sapply(y, weight_func)
  w_dat <- array(sapply(dat, weight_func), dim(dat))
  score <- es_sample(y, dat, w = w_dat)*w_y 
  return (score)
}

################################################################################
# variogram score of order p
#' @rdname scores_sample_multiv_weighted
#' @export
twvs_sample <- function(y, dat, a = -Inf, b = Inf, w_vs = NULL, p = 0.5){
  if (length(a) == 1) {
    a <- rep(a, length(y))
  }else if (length(b) == 1) {
    b <- rep(b, length(y))
  }else if (length(a) != length(y) | length(b) != length(y)) {
    stop("The vector of weight values is not the same length as the vector of realizations")
  }
  v_y <- pmin(pmax(y, a), b)
  v_dat <- sapply(1:nrow(dat), function(m) pmin(pmax(dat[m, ], a), b))
  v_dat <- t(v_dat)
  score <- vs_sample(y = v_y, dat = v_dat, w_vs, p)
  return (score)
}

#' @rdname scores_sample_multiv_weighted
#' @export
owvs_sample <- function(y, dat, a = -Inf, b = Inf, w = NULL, w_vs = NULL, p = 0.5) {
  weight_func <- function(x) as.numeric(x > a & x < b)
  w_y <- sapply(y, weight_func)
  w_dat <- array(sapply(dat, weight_func), dim(dat))
  score <- vs_sample(y, dat, w_vs, p)*w_y # add w argument
  return (score)
}
