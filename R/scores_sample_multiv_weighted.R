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
#' @seealso \link{scores_sample_multiv} for standard (unweighted) scores based on simulated multivariate forecast distributions. \link{scores_sample_univ_weighted} for weighted scores based on simulated univariate forecast distributions
#' 
#' @details
#' In the input matrix \code{dat} each column is expected to represent a sample
#' from the multivariate forecast distribution, the number of rows of \code{dat}
#' thus has to match the length of the observation vector \code{y}, and the
#' number of columns of \code{dat} is the number of simulated samples.
#' 
#' The threshold-weighted scores (\code{\link{twes_sample}}, \code{\link{twmmds_sample}}, 
#' \code{\link{twvs_sample}}) transform \code{y} and \code{dat} using the chaining
#' function \code{chain_func} and then call the relevant unweighted score function
#' (\code{\link{es_sample}}, \code{\link{mmds_sample}}, \code{\link{vs_sample}}). 
#' The outcome-weighted scores (\code{\link{owes_sample}}, \code{\link{owmmds_sample}}, 
#' \code{\link{owvs_sample}}) weight \code{y} and \code{dat} using the weight
#' function \code{weight_func} and then call the relevant unweighted score function
#' (\code{\link{es_sample}}, \code{\link{mmds_sample}}, \code{\link{vs_sample}}). 
#' See the documentation for e.g. \code{\link{es_sample}} for further details.
#'
#' The default weight function used in the weighted scores is 
#' \code{w(z) = 1{a[1] < z[1] < b[1], ..., a[d] < z[d] < b[d]}}, which is equal to one 
#' if \code{z} is in the 'box' defined by the vectors \code{a} and \code{b}, and 
#' is equal to zero otherwise. This weight function emphasises outcomes between 
#' the vectors \code{a} and \code{b}, and is commonly used in practical applications 
#' when interest is on values above a threshold along multiple dimensions.
#'
#' Alternative weight functions can also be employed using the \code{chain_func} 
#' and \code{weight_func} arguments. Computation of the threshold-weighted scores
#' for samples from a predictive distribution requires a chaining function rather 
#' than a weight function. This is why a chaining function is an input for 
#' \code{\link{twes_sample}}, \code{\link{twmmds_sample}}, and \code{\link{twvs_sample}},
#' whereas a weight function is an input for \code{\link{owes_sample}}, 
#' \code{\link{owmmds_sample}}, and \code{\link{owvs_sample}}. 
#'
#' The \code{chain_func} and \code{weight_func} arguments are functions that will 
#' be applied to the elements in \code{y} and \code{dat}. 
#' \code{weight_func} must input a numeric vector of length d, and output a single 
#' numeric value. An error will be returned if \code{weight_func} returns negative values.
#' \code{chain_func} must input a numeric vector of length d, and return a numeric
#' vector of length d.
#'
#' If no custom argument is given for \code{a}, \code{b}, \code{chain_func} or 
#' \code{weight_func}, then all weighted scores are equivalent to the standard 
#' unweighted scores \code{\link{es_sample}}, \code{\link{mmds_sample}}, and
#' \code{\link{vs_sample}}.
#'
#' The \code{w} argument is also present in the unweighted scores.
#' \code{w} is used to weight the draws from the predictive distribution, and does 
#' not weight particular outcomes within the weighted scoring rules. This should not be
#' confused with the \code{weight_func} argument.
#'
#' @return
#' Value of the score. \emph{A lower score indicates a better forecast.}
#' 
#' @references
#' 
#' Allen, S. (2024): 
#' `Weighted scoringRules: Emphasising Particular Outcomes when Evaluating Probabilistic Forecasts', 
#' \emph{Journal of Statistical Software}.
#' \doi{10.18637/jss.v110.i08}
#' 
#' \emph{Threshold-weighted scores}
#' 
#' Allen, S., Ginsbourger, D. and J. Ziegel (2023): 
#' `Evaluating forecasts for high-impact events using transformed kernel scores', 
#' \emph{SIAM/ASA Journal on Uncertainty Quantification} 11, 906-940.
#' \doi{10.1137/22M1532184}
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
#' @examples
#' \dontrun{
#' d <- 3  # number of dimensions
#' m <- 10  # number of samples from multivariate forecast distribution
#' 
#' # parameters for multivariate normal example
#' mu0 <- rep(0, d)
#' mu <- rep(1, d)
#' S0 <- S <- diag(d)
#' S0[S0==0] <- 0.2
#' S[S==0] <- 0.1
#' 
#' # generate samples from multivariate normal distributions
#' obs <- drop(mu0 + rnorm(d) %*% chol(S0))
#' sample_fc <- replicate(m, drop(mu + rnorm(d) %*% chol(S)))
#' 
#' # if no additional parameters are provided, the weighted scores are the same as
#' # the unweighted scores:
#' es_sample(y = obs, dat = sample_fc) # energy score
#' twes_sample(y = obs, dat = sample_fc) # threshold-weighted energy score
#' owes_sample(y = obs, dat = sample_fc) # outcome-weighted energy score
#' 
#' mmds_sample(y = obs, dat = sample_fc) # Gaussian kernel score
#' twmmds_sample(y = obs, dat = sample_fc) # threshold-weighted Gaussian kernel score
#' owmmds_sample(y = obs, dat = sample_fc) # outcome-weighted Gaussian kernel score
#' 
#' vs_sample(y = obs, dat = sample_fc) # variogram score
#' twvs_sample(y = obs, dat = sample_fc) # threshold-weighted variogram score
#' owvs_sample(y = obs, dat = sample_fc) # outcome-weighted variogram score
#' 
#'
#' # the outcome-weighted scores are undefined if none of dat are between a and b
#' # this can lead to NaNs in some of the scores calculated below, particularly
#' # if the thresholds are extreme, or if the dimension is large
#'
#'
#' # emphasise outcomes greater than 0 in all dimensions
#' twes_sample(y = obs, dat = sample_fc, a = 0)
#' owes_sample(y = obs, dat = sample_fc, a = 0)
#' twmmds_sample(y = obs, dat = sample_fc, a = 0)
#' owmmds_sample(y = obs, dat = sample_fc, a = 0)
#' twvs_sample(y = obs, dat = sample_fc, a = 0)
#' owvs_sample(y = obs, dat = sample_fc, a = 0)
#'
#' # this can also be done more explicitly by setting a = rep(0, d)
#' twes_sample(y = obs, dat = sample_fc, a = rep(0, d))
#' owes_sample(y = obs, dat = sample_fc, a = rep(0, d))
#'
#' # a should also be specified fully if the threshold changes in each dimension
#' a <- rnorm(d)
#' twes_sample(y = obs, dat = sample_fc, a = a)
#' owes_sample(y = obs, dat = sample_fc, a = a)
#' twmmds_sample(y = obs, dat = sample_fc, a = a)
#' owmmds_sample(y = obs, dat = sample_fc, a = a)
#' twvs_sample(y = obs, dat = sample_fc, a = a)
#' owvs_sample(y = obs, dat = sample_fc, a = a)
#'
#' # emphasise outcomes smaller than 0 in all dimensions
#' twes_sample(y = obs, dat = sample_fc, b = 0)
#' owes_sample(y = obs, dat = sample_fc, b = 0)
#' twmmds_sample(y = obs, dat = sample_fc, b = 0)
#' owmmds_sample(y = obs, dat = sample_fc, b = 0)
#' twvs_sample(y = obs, dat = sample_fc, b = 0)
#' owvs_sample(y = obs, dat = sample_fc, b = 0)
#'
#' # emphasise outcomes between (-1, -1, -1) and (1, 1, 1)
#' twes_sample(y = obs, dat = sample_fc, a = -1, b = 1)
#' owes_sample(y = obs, dat = sample_fc, a = -1, b = 1)
#' twmmds_sample(y = obs, dat = sample_fc, a = -1, b = 1)
#' owmmds_sample(y = obs, dat = sample_fc, a = -1, b = 1)
#' twvs_sample(y = obs, dat = sample_fc, a = -1, b = 1)
#' owvs_sample(y = obs, dat = sample_fc, a = -1, b = 1)
#'
#' # emphasise outcomes between (-2, 0, -1) and (0, 2, 1)
#' a <- c(-2, 0, -1)
#' b <- c(0, 2, 1)
#' twes_sample(y = obs, dat = sample_fc, a = a, b = b)
#' owes_sample(y = obs, dat = sample_fc, a = a, b = b)
#' twmmds_sample(y = obs, dat = sample_fc, a = a, b = b)
#' owmmds_sample(y = obs, dat = sample_fc, a = a, b = b)
#' twvs_sample(y = obs, dat = sample_fc, a = a, b = b)
#' owvs_sample(y = obs, dat = sample_fc, a = a, b = b)
#'
#'
#' # values of a cannot be larger than the corresponding values of b
#' twes_sample(y = obs, dat = sample_fc, a = c(0, 0, 0), b = c(0, 0, 1))
#' twes_sample(y = obs, dat = sample_fc, a = c(0, 0, 0), b = c(0, 0, 0)) # error
#' twes_sample(y = obs, dat = sample_fc, a = c(0, 0, 0), b = c(1, 1, -1)) # error
#' 
#' # a and b must be of the same length (and of the same length as y)
#' owmmds_sample(y = obs, dat = sample_fc, a = c(0, 0), b = 1) # error
#' owmmds_sample(y = obs, dat = sample_fc, a = c(0, 0), b = c(1, 1)) # error
#'
#'
#' # alternative custom weight and chaining functions can also be used
#'
#' # Example 1: the default weight function with an alternative chaining function
#' # the default weight function is 
#' # w(z) = 1{a[1] < z[1] < b[1], ..., a[d] < z[d] < b[d]}
#' # the default chaining function is 
#' # v(z) = (min(max(z[1], a[1]), b[1]), ..., min(max(z[d], a[d]), b[d]))
#' a <- -2
#' b <- 2
#' weight_func <- function(x) as.numeric(all(x > a & x < b))
#' chain_func <- function(x) pmin(pmax(x, a), b)
#' owes_sample(y = obs, dat = sample_fc, a = a, b = b)
#' owes_sample(y = obs, dat = sample_fc, weight_func = weight_func)
#' twes_sample(y = obs, dat = sample_fc, a = a, b = b)
#' twes_sample(y = obs, dat = sample_fc, chain_func = chain_func)
#'
#' # consider an alternative chaining function: v(z) = z if w(z) = 1, else v(z) = 0
#' chain_func <- function(x) x*weight_func(x)
#' twes_sample(y = obs, dat = sample_fc, chain_func = chain_func)
#'
#'
#' # Example 2: a mulivariate Gaussian weight function with mean vector mu and 
#' # diagonal covariance matrix sigma
#' mu <- rep(0, d)
#' sigma <- diag(d)
#' weight_func <- function(x) prod(pnorm(x, mu, diag(sigma)))
#' # the corresponding chaining function is
#' chain_func <- function(x){
#'  (x - mu)*pnorm(x, mu, diag(sigma)) + (diag(sigma)^2)*dnorm(x, mu, diag(sigma))
#' }
#'
#' owvs_sample(y = obs, dat = sample_fc, a = mu)
#' owvs_sample(y = obs, dat = sample_fc, weight_func = weight_func)
#' twvs_sample(y = obs, dat = sample_fc, a = mu)
#' twvs_sample(y = obs, dat = sample_fc, chain_func = chain_func)
#'}
#' 
#' @name scores_sample_multiv_weighted
NULL

################################################################################
# threshold-weighted energy score
#' @rdname scores_sample_multiv_weighted
#' @export
twes_sample <- function(y, dat, a = -Inf, b = Inf, chain_func = NULL , w = NULL) {
  if (identical(length(a), 1L)) {
    a <- rep(a, length(y))
  }
  if (identical(length(b), 1L)) {
    b <- rep(b, length(y))
  }
  
  input <- list(lower = a, upper = b, v = chain_func, y = y)
  check_mv_weight(input)
  
  if (is.null(chain_func)) chain_func <- function(x) pmin(pmax(x, a), b)
  
  v_y <- chain_func(y)
  v_dat <- apply(dat, 2, chain_func)
  score <- es_sample(y = v_y, dat = v_dat, w)
  return(score)
}

################################################################################
# outcome-weighted energy score
#' @rdname scores_sample_multiv_weighted
#' @export
owes_sample <- function(y, dat, a = -Inf, b = Inf, weight_func = NULL, w = NULL) {
  if (identical(length(a), 1L)) {
    a <- rep(a, length(y))
  }
  if (identical(length(b), 1L)) {
    b <- rep(b, length(y))
  }
  
  input <- list(lower = a, upper = b, w = weight_func, y = y)
  check_mv_weight(input)
  
  if (is.null(weight_func)) weight_func <- function(x) as.numeric(all(x > a & x < b))
  
  w_y <- weight_func(y)
  w_dat <- apply(dat, 2, weight_func)
  if (is.null(w)) {
    w <- w_dat
  }else {
    w <- w*w_dat
  }
  if (identical(sum(w), 0)) {
    score <- NaN
  } else { 
    w <- w/sum(w)
    score <- es_sample(y, dat, w = w)*w_y 
  }
  return(score)
}

################################################################################
# threshold-weighted MMD score
#' @rdname scores_sample_multiv_weighted
#' @export
twmmds_sample <- function(y, dat, a = -Inf, b = Inf, chain_func = NULL, w = NULL) {
  if (identical(length(a), 1L)) {
    a <- rep(a, length(y))
  }
  if (identical(length(b), 1L)) {
    b <- rep(b, length(y))
  }
  
  input <- list(lower = a, upper = b, v = chain_func, y = y)
  check_mv_weight(input)
  
  if (is.null(chain_func)) chain_func <- function(x) pmin(pmax(x, a), b)
  
  v_y <- chain_func(y)
  v_dat <- apply(dat, 2, chain_func)
  score <- mmds_sample(y = v_y, dat = v_dat, w)
  return(score)
}

################################################################################
# outcome-weighted MMD score
#' @rdname scores_sample_multiv_weighted
#' @export
owmmds_sample <- function(y, dat, a = -Inf, b = Inf, weight_func = NULL, w = NULL) {
  if (identical(length(a), 1L)) {
    a <- rep(a, length(y))
  }
  if (identical(length(b), 1L)) {
    b <- rep(b, length(y))
  }
  
  input <- list(lower = a, upper = b, w = weight_func, y = y)
  check_mv_weight(input)
  
  if (is.null(weight_func)) weight_func <- function(x) as.numeric(all(x > a & x < b))
  
  w_y <- weight_func(y)
  w_dat <- apply(dat, 2, weight_func)
  if (is.null(w)) {
    w <- w_dat
  }else {
    w <- w*w_dat
  }
  if (identical(sum(w), 0)) {
    score <- NaN
  } else { 
    w <- w/sum(w)
    score <- mmds_sample(y, dat, w = w)*w_y 
  }
  return(score)
}

################################################################################
# threshold-weighted variogram score of order p
#' @rdname scores_sample_multiv_weighted
#' @export
twvs_sample <- function(y, dat, a = -Inf, b = Inf, chain_func = NULL, w = NULL, w_vs = NULL, p = 0.5){
  if (identical(length(a), 1L)) {
    a <- rep(a, length(y))
  }
  if (identical(length(b), 1L)) {
    b <- rep(b, length(y))
  }
  
  input <- list(lower = a, upper = b, v = chain_func, y = y)
  check_mv_weight(input)
  
  if (is.null(chain_func)) chain_func <- function(x) pmin(pmax(x, a), b)
  
  v_y <- chain_func(y)
  v_dat <- apply(dat, 2, chain_func)
  score <- vs_sample(y = v_y, dat = v_dat, w = w, w_vs = w_vs, p = p)
  return(score)
}

################################################################################
# outcome-weighted variogram score of order p
#' @rdname scores_sample_multiv_weighted
#' @export
owvs_sample <- function(y, dat, a = -Inf, b = Inf, weight_func = NULL, w = NULL, w_vs = NULL, p = 0.5) {
  if (identical(length(a), 1L)) {
    a <- rep(a, length(y))
  }
  if (identical(length(b), 1L)) {
    b <- rep(b, length(y))
  }
  
  input <- list(lower = a, upper = b, w = weight_func, y = y)
  check_mv_weight(input)
  
  if (is.null(weight_func)) weight_func <- function(x) as.numeric(all(x > a & x < b))
  
  w_y <- weight_func(y)
  w_dat <- apply(dat, 2, weight_func)
  if (is.null(w)) {
    w <- w_dat
  }else {
    w <- w*w_dat
  }
  if (identical(sum(w), 0)) {
    score <- NaN
  } else { 
    w <- w/sum(w)
    score <- vs_sample(y, dat, w = w, w_vs = w_vs, p = p)*w_y 
  }
  return(score)
}

################################################################################
# checks for the weight and chaining functions in the multivariate weighted scoring rules
check_mv_weight <- function(input) {
  a <- input$lower
  b <- input$upper
  w <- input$w
  v <- input$v
  y <- input$y
  
  if (!identical(length(a), length(b))) {
    stop("a and b do not have the same length.")
  } else if (!identical(length(a), length(y))) {
    stop("a and b do not have the same length as the realized values.")
  } else if (!is.numeric(a) || !is.numeric(b)) {
    stop("The lower and upper bounds in the weight function must be numeric.")
  } else if (any(a > b) || all(a == b)) {
    stop("The lower bound in the weight function is not smaller than the upper bound.")
  }
  
  if (!is.null(w)) {
    if (!(all(is.infinite(a)) && all(is.infinite(b)))) {
      warning("The arguments 'a' and/or 'b' have been given in addition to a custom weight function. 'a' and 'b' will be ignored.")
    }
    if (!is.function(w)) {
      stop("The weight function must be of type 'function'.")
    } else {
      w_y <- w(y)
    }
    if (!is.numeric(w_y)) {
      stop("The weight function does not return a single numerical value.")
    }else if (!identical(length(w_y), 1L)) {
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
    }else if (!identical(length(v_y), length(y))) {
      stop("The chaining function does not return a numeric vector of the same length as y.")
    }
  }
}

