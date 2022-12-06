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
twes_sample <- function(y, dat, v = identity, w = NULL){

  v_y <- v(y)
  check.chaining(y, v_y)
  v_dat <- apply(dat, 2, v)
  check.chaining(dat, v_dat)
  
  out <- es_sample(y = v_y, dat = v_dat, w = w)
  
  return(out)
}

################################################################################
# MMD score
#' @rdname scores_sample_multiv_weighted
#' @export
twmmds_sample <- function(y, dat, v = identity, w = NULL){
  
  v_y <- v(y)
  check.chaining(y, v_y)
  v_dat <- apply(dat, 2, v)
  check.chaining(dat, v_dat)
  
  out <- es_sample(y = v_y, dat = v_dat, w = w)
  
  return(out)
}

################################################################################
# variogram score of order p
#' @rdname scores_sample_multiv_weighted
#' @export
twvs_sample <- function(y, dat, v = identity, w_vs = NULL, p = 0.5){
  # y: realised values (numeric vector of length d)
  # dat: numeric matrix of data (columns are simulation draws from multivariate forecast distribution)
  # v: chaining function to transform the forecasts and observations (function)
  # w_vs: numeric matrix of weights for dat used in the variogram score. This matrix must be square and symmetric, with all elements being non-negative.
  # p: order of variogram score. Standard choices include p = 1 and p = 0.5.
  
  v_y <- v(y)
  check.chaining(y, v_y)
  v_dat <- apply(dat, 2, v)
  check.chaining(dat, v_dat)
  
  out <- vs_sample(y = v_y, dat = v_dat, w_vs = w_vs, p = p)
  
  return(out)
}

