#' Weighted Multivariate Scoring Rules for Simulated Forecast Distributions
#' 
#' Compute weighted versions of multivariate scores \eqn{S(y, dat)}, where \eqn{S} is a
#' proper scoring rule, \eqn{y} is a d-dimensional realization vector and 
#' \eqn{dat} is a simulated sample of multivariate forecasts. The weighted scores allow 
#' particular outcomes of interest to be emphasised during forecast evaluation.
#' Threshold-weighted and outcome-weighted versions of three multivariate scores are 
#' available: the energy score, a score based on a Gaussian kernel (\link{mmds_sample}, 
#' see details below) and the variogram score of order \eqn{p}.
#' 
#' @param y realized values (numeric vector of length d).
#' @param dat numeric matrix of data
#' (columns are simulation draws from multivariate forecast distribution).
#' @param a numeric vector of of length d containing lower bounds for the indicator 
#' weight function \code{w(z) = 1{a[1] < z[1] < b[1], ..., a[d] < z[d] < b[d]}}.
#' @param b numeric vector of of length d containing upper bounds for the indicator 
#' weight function \code{w(z) = 1{a[1] < z[1] < b[1], ..., a[d] < z[d] < b[d]}}.
#' @param chain_func function used to target particular outcomes in the threshold-weighted scores; 
#' the default corresponds to the weight function above.
#' @param weight_func function used to target particular outcomes in the outcome-weighted scores; 
#' the default corresponds to the weight function above.
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
# threshold-weighted energy score
#' @rdname scores_sample_multiv_weighted
#' @export
twes_sample <- function(y, dat, a = -Inf, b = Inf, chain_func = NULL, w = NULL) {
  if (length(a) == 1) {
    a <- rep(a, length(y))
  }
  if (length(b) == 1) {
    b <- rep(b, length(y))
  }
  input <- list(lower = a, upper = b, v = chain_func, y = y)
  check_mv_weight(input)
  if (is.null(chain_func)) {
    chain_func <- function(x) pmin(pmax(x, a), b) 
  }
  v_y <- chain_func(y)
  v_dat <- sapply(1:ncol(dat), function(m) chain_func(dat[, m]))
  score <- es_sample(y = v_y, dat = v_dat, w)
  return(score)
}

################################################################################
# outcome-weighted energy score
#' @rdname scores_sample_multiv_weighted
#' @export
owes_sample <- function(y, dat, a = -Inf, b = Inf, weight_func = NULL, w = NULL) {
  if (length(a) == 1) {
    a <- rep(a, length(y))
  }
  if (length(b) == 1) {
    b <- rep(b, length(y))
  }
  input <- list(lower = a, upper = b, w = weight_func, y = y)
  check_mv_weight(input)
  if (is.null(weight_func)) {
    weight_func <- function(x) as.numeric(all(x > a & x < b))
  }
  w_y <- weight_func(y)
  w_dat <- sapply(1:ncol(dat), function(m) weight_func(dat[, m]))
  if (is.null(w)) {
    w <- w_dat
  }else {
    w <- w*w_dat
  }
  score <- es_sample(y, dat, w = w)*w_y 
  return(score)
}

################################################################################
# threshold-weighted MMD score
#' @rdname scores_sample_multiv_weighted
#' @export
twmmds_sample <- function(y, dat, a = -Inf, b = Inf, chain_func = NULL, w = NULL) {
  if (length(a) == 1) {
    a <- rep(a, length(y))
  }
  if (length(b) == 1) {
    b <- rep(b, length(y))
  }
  input <- list(lower = a, upper = b, v = chain_func, y = y)
  check_mv_weight(input)
  if (is.null(chain_func)) {
    chain_func <- function(x) pmin(pmax(x, a), b) 
  }
  v_y <- chain_func(y)
  v_dat <- sapply(1:ncol(dat), function(m) chain_func(dat[, m]))
  score <- mmds_sample(y = v_y, dat = v_dat, w)
  return(score)
}

################################################################################
# outcome-weighted MMD score
#' @rdname scores_sample_multiv_weighted
#' @export
owmmds_sample <- function(y, dat, a = -Inf, b = Inf, weight_func = NULL, w = NULL) {
  if (length(a) == 1) {
    a <- rep(a, length(y))
  }
  if (length(b) == 1) {
    b <- rep(b, length(y))
  }
  input <- list(lower = a, upper = b, w = weight_func, y = y)
  check_mv_weight(input)
  if (is.null(weight_func)) {
    weight_func <- function(x) as.numeric(all(x > a & x < b))
  }
  w_y <- weight_func(y)
  w_dat <- sapply(1:ncol(dat), function(m) weight_func(dat[, m]))
  if (is.null(w)) {
    w <- w_dat
  }else {
    w <- w*w_dat
  }
  score <- mmds_sample(y, dat, w = w)*w_y 
  return(score)
}

################################################################################
# threshold-weighted variogram score of order p
#' @rdname scores_sample_multiv_weighted
#' @export
twvs_sample <- function(y, dat, a = -Inf, b = Inf, chain_func = NULL, w = NULL, w_vs = NULL, p = 0.5){
  if (length(a) == 1) {
    a <- rep(a, length(y))
  }
  if (length(b) == 1) {
    b <- rep(b, length(y))
  }
  input <- list(lower = a, upper = b, v = chain_func, y = y)
  check_mv_weight(input)
  if (is.null(chain_func)) {
    chain_func <- function(x) pmin(pmax(x, a), b) 
  }
  v_y <- chain_func(y)
  v_dat <- sapply(1:ncol(dat), function(m) chain_func(dat[, m]))
  score <- vs_sample(y = v_y, dat = v_dat, w = w, w_vs = w_vs, p = p)
  return(score)
}

################################################################################
# outcome-weighted variogram score of order p
#' @rdname scores_sample_multiv_weighted
#' @export
owvs_sample <- function(y, dat, a = -Inf, b = Inf, weight_func = NULL, w = NULL, w_vs = NULL, p = 0.5) {
  if (length(a) == 1) {
    a <- rep(a, length(y))
  }
  if (length(b) == 1) {
    b <- rep(b, length(y))
  }
  input <- list(lower = a, upper = b, w = weight_func, y = y)
  check_mv_weight(input)
  if (is.null(weight_func)) {
    weight_func <- function(x) as.numeric(all(x > a & x < b))
  }
  w_y <- weight_func(y)
  w_dat <- sapply(1:ncol(dat), function(m) weight_func(dat[, m]))
  if (is.null(w)) {
    w <- w_dat
  }else {
    w <- w*w_dat
  }
  score <- vs_sample(y, dat, w = w, w_vs = w_vs, p = p)*w_y 
  return(score)
}

################################################################################
# checks for the weight and chaining functions in the multivariate weighted scoring rules
check_mv_weight <- function(input) {
  a <- input$lower
  b <- input$upper
  w <- input$w
  v <- input$v
  
  if (length(a) != length(b)) {
    stop("a and b do not have the same length.")
  } else if (length(a) != length(y)) {
    stop("a and b do not have the same length as the realized values.")
  } else if (!is.numeric(a) | !is.numeric(b)) {
    stop("The lower and upper bounds in the weight function must be numeric.")
  } else if (any(a > b) | all(a == b)) {
    stop("The lower bound in the weight function is not smaller than the upper bound.")
  }
  
  if (!is.null(w)) {
    if (!is.function(w)) {
      stop("The weight function must be of type 'function'.")
    } else {
      w_y <- w(y)
    }
    if (!is.numeric(w_y)) {
      stop("The weight function does not return a single numerical value.")
    }else if (length(w_y) != 1) {
      stop("The weight function does not return a single numerical value.")
    } else if (w_y < 0) {
      stop("The weight function returns negative weights.")
    }
  }
  
  if (!is.null(v)) {
    if (!is.function(v)) {
      stop("The chaining function must be of type 'function'.")
    } else {
      v_y <- v(y)
    }
    if (any(!is.numeric(v_y))) {
      stop("The chaining function does not return a numeric vector.")
    }else if (length(v_y) != length(y)) {
      stop("The chaining function does not return a numeric vector of the same length as y.")
    }
  }
}

