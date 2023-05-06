#' Weighted Scoring Rules for Simulated Forecast Distributions (experimental)
#' 
#' Calculate weighted scores given observations and draws from univariate predictive distributions.
#' The weighted scoring rules that are available are the threshold-weighted CRPS, outcome-weighted CRPS, 
#' and conditional and censored likelihood scores. Note that the functions 
#' documented here are a new experimental feature of the package, and feedback is highly welcome.
#' 
#' @param y vector of realized values.
#' @param dat vector or matrix (depending on \code{y}; see details)
#'  of simulation draws from forecast distribution. 
#' @param a numeric lower bound for the indicator weight function \code{w(z) = 1{a < z < b}}.
#' @param b numeric upper bound for the indicator weight function \code{w(z) = 1{a < z < b}}.
#' @param chain_func function used to target particular outcomes in the threshold-weighted CRPS; 
#' the default corresponds to the weight function \code{w(z) = 1{a < z < b}}.
#' @param weight_func function used to target particular outcomes in the outcome-weighted CRPS; 
#' the default corresponds to the weight function \code{w(z) = 1{a < z < b}}.
#' @param w optional; vector or matrix (matching \code{dat}) of ensemble weights. 
#'  Note that these weights are not used in the weighted scoring rules; see details.
#' @param bw optional; vector (matching \code{y}) of bandwidths for kernel density
#' estimation for \code{\link{clogs_sample}}; see details.
#' @param show_messages logical; display of messages (does not affect
#'  warnings and errors).
#' @param cens logical; if TRUE, \code{\link{clogs_sample}} returns the censored
#'  likelihood score; if FALSE, \code{\link{clogs_sample}} returns the conditional
#'  likelihood score.
#'  
#' @return
#' Value of the score. \emph{A lower score indicates a better forecast.}
#' 
#' @references
#' \emph{Threshold-weighted CRPS:}
#' 
#' Gneiting, T. and R. Ranjan (2011): 
#' `Comparing density forecasts using threshold-and quantile-weighted scoring rules', 
#' \emph{Journal of Business & Economic Statistics} 29, 411-422. 
#' \doi{10.1198/jbes.2010.08110}
#' 
#' Allen, S., Ginsbourger, D. and J. Ziegel (2022): 
#' `Evaluating forecasts for high-impact events using transformed kernel scores', 
#' \emph{arXiv preprint} arXiv:2202.12732.
#' \doi{10.48550/arXiv.2202.12732}
#'  
#' \emph{Outcome-weighted CRPS:}
#'  
#' Holzmann, H. and B. Klar (2017):
#' `Focusing on regions of interest in forecast evaluation',
#' \emph{Annals of Applied Statistics} 11, 2404-2431. 
#' \doi{10.1214/17-AOAS1088}
#' 
#' \emph{Conditional and censored likelihood scores:}
#' 
#' Diks, C., Panchenko, V. and D. Van Dijk (2011):
#' `Likelihood-based scoring rules for comparing density forecasts in tails',
#' \emph{Journal of Econometrics} 163, 215-230.
#' \doi{10.1016/j.jeconom.2011.04.001}
#'  
#' @author Sam Allen
#' 
#' @seealso \link{scores_sample_univ} for standard (unweighted) scores based on simulated forecast distributions. \link{scores_sample_multiv_weighted} for weighted scores based on simulated multivariate forecast distributions.
#' 
#' @details 
#' For a vector \code{y} of length n, \code{dat} should be given as a matrix
#' with n rows. If \code{y} has length 1, then \code{dat} may be a vector.
#' 
#' \code{\link{twcrps_sample}} transforms \code{y} and \code{dat} using the chaining
#' function \code{chain_func} and then calls \code{\link{crps_sample}}. 
#' \code{\link{owcrps_sample}} weights \code{y} and \code{dat} using the weight function
#' \code{weight_func} and then calls \code{\link{crps_sample}}. 
#' See the documentation for \code{\link{crps_sample}} for further details.
#' 
#' The default weight function used in the weighted scores is \code{w(z) = 1{a < z < b}}, 
#' which is equal to one if \code{z} is between \code{a} and \code{b}, and zero otherwise.
#' This weight function emphasises outcomes between \code{a} and \code{b}, and is 
#' commonly used in practical applications when interest is on values above a threshold
#' (set \code{b = Inf} and \code{a} equal to the threshold) or below a threshold 
#' (set \code{a = -Inf} and \code{b} equal to the threshold). 
#' 
#' Alternative weight functions can also be employed using the \code{chain_func} 
#' and \code{weight_func} arguments to \code{\link{twcrps_sample}} and \code{\link{owcrps_sample}},
#' respectively. Computation of the threshold-weighted CRPS for samples from a predictive distribution 
#' requires a chaining function rather than a weight function. This is why a chaining 
#' function is an input for \code{\link{twcrps_sample}} whereas a weight function is an 
#' input for \code{\link{owcrps_sample}}. Since \code{\link{clogs_sample}} requires 
#' kernel density estimation to approximate the forecast density, it cannot readily
#' be calculated for arbitrary weight functions, and is thus only available for 
#' the canonical weight function \code{w(z) = 1{a < z < b}}.
#' 
#' The \code{chain_func} and \code{weight_func} arguments are functions that will 
#' be applied to the vector \code{y} and the columns of \code{dat}. It is assumed
#' that these functions are vectorised. Both functions must take a vector as an input
#' and output a vector of the same length, containing the weight (for \code{weight_func}) 
#' or transformed value (for \code{chain_func}) corresponding to each element in the 
#' input vector. An error will be returned if \code{weight_func} returns
#' negative values, and a warning message will appear if \code{chain_func} is 
#' not increasing. 
#' 
#' If no custom argument is given for \code{a}, \code{b}, \code{chain_func} or 
#' \code{weight_func}, then both \code{\link{twcrps_sample}} and \code{\link{owcrps_sample}} 
#' are equivalent to the standard unweighted \code{\link{crps_sample}}, and 
#' \code{\link{clogs_sample}} is equivalent to \code{\link{logs_sample}}. 
#' 
#' The \code{w} argument is also present in the unweighted scores (e.g. \code{\link{crps_sample}}).
#' \code{w} is used to weight the draws from the predictive distribution, and does 
#' not weight particular outcomes within the weighted scoring rules. This should not be
#' confused with the \code{weight_func} argument, which is used within the weighted scores.
#' 
#' @examples
#' \dontrun{
#' 
#' y <- rnorm(10)
#' sample_fc <- matrix(rnorm(100), nrow = 10)
#' 
#' crps_sample(y = y, dat = sample_fc)
#' twcrps_sample(y = y, dat = sample_fc)
#' owcrps_sample(y = y, dat = sample_fc)
#' 
#' logs_sample(y = y, dat = sample_fc)
#' clogs_sample(y = y, dat = sample_fc)
#' clogs_sample(y = y, dat = sample_fc, cens = FALSE)
#' 
#' # emphasise outcomes above 0
#' twcrps_sample(y = y, dat = sample_fc, a = 0)
#' owcrps_sample(y = y, dat = sample_fc, a = 0)
#' clogs_sample(y = y, dat = sample_fc, a = 0)
#' clogs_sample(y = y, dat = sample_fc, a = 0, cens = FALSE)
#' 
#' # emphasise outcomes below 0
#' twcrps_sample(y = y, dat = sample_fc, b = 0)
#' owcrps_sample(y = y, dat = sample_fc, b = 0)
#' clogs_sample(y = y, dat = sample_fc, b = 0) 
#' 
#' # emphasise outcomes between -1 and 1
#' twcrps_sample(y = y, dat = sample_fc, a = -1, b = 1)
#' owcrps_sample(y = y, dat = sample_fc, a = -1, b = 1)
#' clogs_sample(y = y, dat = sample_fc, a = -1, b = 1)
#' 
#' 
#' # a must be smaller than b 
#' twcrps_sample(y = y, dat = sample_fc, a = 1, b = -1) # error
#' owcrps_sample(y = y, dat = sample_fc, a = 0, b = 0) # error
#' clogs_sample(y = y, dat = sample_fc, a = 10, b = 9) # error
#' 
#' # a and b must be single numeric values (not vectors)
#' twcrps_sample(y = y, dat = sample_fc, a = rnorm(10)) # error
#' 
#' 
#' # the owCRPS is not well-defined if none of dat are between a and b
#' y <- rnorm(10)
#' sample_fc <- matrix(runif(100, -5, 1), nrow = 10)
#' owcrps_sample(y = y, dat = sample_fc, a = 1)
#' # the twCRPS is zero if none of y and dat are between a and b
#' twcrps_sample(y = y, dat = sample_fc, a = 1) 
#' 
#' 
#' # alternative custom weight and chaining functions can also be used
#' 
#' # Example 1: a Gaussian weight function with location mu and scale sigma
#' mu <- 0
#' sigma <- 0.5
#' weight_func <- function(x) pnorm(x, mu, sigma)
#' # a corresponding chaining function is
#' chain_func <- function(x) (x - mu)*pnorm(x, mu, sigma) + (sigma^2)*dnorm(x, mu, sigma)
#' 
#' x <- seq(-2, 2, 0.01)
#' plot(x, weight_func(x), type = "l") # positive outcomes are given higher weight
#' plot(x, chain_func(x), type = "l") 
#' 
#' owcrps_sample(y = y, dat = sample_fc, a = mu)
#' owcrps_sample(y = y, dat = sample_fc, weight_func = weight_func)
#' twcrps_sample(y = y, dat = sample_fc, a = mu)
#' twcrps_sample(y = y, dat = sample_fc, chain_func = chain_func)
#' 
#' 
#' # Example 2: a sigmoid (or logistic) weight function with location mu and scale sigma
#' weight_func <- function(x) plogis(x, mu, sigma)
#' chain_func <- function(x) sigma*log(exp((x - mu)/sigma) + 1)
#' 
#' x <- seq(-2, 2, 0.01)
#' plot(x, weight_func(x), type = "l") # positive outcomes are given higher weight
#' plot(x, chain_func(x), type = "l") 
#' 
#' owcrps_sample(y = y, dat = sample_fc, a = mu)
#' owcrps_sample(y = y, dat = sample_fc, weight_func = weight_func)
#' twcrps_sample(y = y, dat = sample_fc, a = mu)
#' twcrps_sample(y = y, dat = sample_fc, chain_func = chain_func)
#' 
#' 
#' # Example 3: the weight function w(z) = 1{z < a or z > b}
#' a <- -1
#' b <- 1
#' weight_func <- function(x) as.numeric(x < a | x > b)
#' chain_func <- function(x) (x < a)*(x - a) + (x > b)*(x - b) + a
#' 
#' x <- seq(-2, 2, 0.01)
#' plot(x, weight_func(x), type = "l")
#' plot(x, chain_func(x), type = "l")
#' 
#' owcrps_sample(y = y, dat = sample_fc, weight_func = weight_func)
#' twcrps_sample(y = y, dat = sample_fc, chain_func = chain_func)
#' twcrps_sample(y = y, dat = sample_fc, b = -1) + twcrps_sample(y = y, dat = sample_fc, a = 1)
#' crps_sample(y = y, dat = sample_fc) - twcrps_sample(y = y, dat = sample_fc, a = -1, b = 1)
#' }
#' 
#' @name scores_sample_univ_weighted
NULL

################################################################################
# threshold-weighted CRPS

#' @rdname scores_sample_univ_weighted
#' @export
twcrps_sample <- function (y, dat, a = -Inf, b = Inf, chain_func = function(x) pmin(pmax(x, a), b),
                           w = NULL, show_messages = TRUE) {
  input <- list(lower = a, upper = b, v = chain_func, y = y)
  check_weight(input)
  v_y <- chain_func(y)
  if (is.vector(dat)) {
    v_dat <- sapply(dat, chain_func)
  } else {
    v_dat <- apply(dat, 2, chain_func)
  }
  crps_sample(y = v_y, dat = v_dat, method = "edf", w = w, show_messages = show_messages)
}

################################################################################
# outcome-weighted CRPS
#' @rdname scores_sample_univ_weighted
#' @export
owcrps_sample <- function (y, dat, a = -Inf, b = Inf, weight_func = function(x) as.numeric(x > a & x < b),
                           w = NULL, show_messages = TRUE) {
  input <- list(lower = a, upper = b, w = weight_func, y = y)
  check_weight(input)
  w_y <- weight_func(y)
  if (is.vector(dat)) {
    w_dat <- sapply(dat, weight_func)
  } else {
    w_dat <- apply(dat, 2, weight_func)
  }
  if (is.null(w)) {
    w <- w_dat
  }else {
    if (identical(dim(w), dim(dat))) {
      w <- w*w_dat
    }else {
      stop("The dimensions of w do not match the dimensions of dat.")
    }
  }
  crps_sample(y, dat, method = "edf", w = w, show_messages = show_messages)*w_y
}

################################################################################
# conditional and censored likelihood scores
#' @rdname scores_sample_univ_weighted
#' @export
clogs_sample <- function (y, dat, a = -Inf, b = Inf, bw = NULL, show_messages = FALSE, cens = TRUE) {
  input <- list(lower = a, upper = b)
  check_weight(input)
  input <- list(y = y, dat = dat)
  input$bw <- bw
  if (show_messages)
    message("Using the log score with kernel density estimation tends to be fragile -- see KLTG (2021) for details.")
  if (identical(length(y), 1L) && is.vector(dat)) {
    check_sample(input)
    if (is.null(bw)) bw <- bw.nrd(dat)
    n <- length(dat)
    w <- rep(1 / n, n)
    s <- rep(bw, n)
    int_wf <- pmixnC(dat, s, b) - pmixnC(dat, s, a)
    ls <- lsmixnC(w, dat, s, y)
  } else {
    check_sample2(input)
    if (is.null(bw)) bw <- apply(dat, 1, bw.nrd)
    w <- rep(1, ncol(dat))
    s <- matrix(bw, nrow(dat), ncol(dat))
    int_wf <- sapply(seq_along(y), function(i) pmixnC(dat[i, ], s[i, ], b) - pmixnC(dat[i, ], s[i, ], a))
    ls <- sapply(seq_along(y), function(i) lsmixnC(w, dat[i, ], s[i, ], y[i]))
  }
  w_y <- as.numeric(y > a & y < b)
  score <- w_y*ls + log(int_wf^w_y)
  if (cens) {
    score <- score - log(int_wf^w_y) - log((1 - int_wf)^(1 - w_y))
  }
  return (score)
}

################################################################################
# checks for the weight and chaining functions in the weighted scoring rules
check_weight <- function(input, tol = 1e-15) {
  a <- input$lower
  b <- input$upper
  w <- input$w
  v <- input$v
  y <- sort(input$y)
  
  if (!is.numeric(a) || !is.numeric(b)) {
    stop("The lower and upper bounds in the weight function must be single numeric values.")
  } else if (length(a) > 1 || length(b) > 1) {
    stop("The lower and upper bounds in the weight function must be single numeric values.")
  } else if (a >= b) {
    stop("The lower bound in the weight function is not smaller than the upper bound.")
  }
  
  if (!is.null(w)) {
    if (!is.function(w)) {
      stop("The weight function must be of type 'function'.")
    } else {
      w_y <- w(y)
    }
    if (any(!is.numeric(w_y))) {
      stop("The weight function does not return a numeric vector.")
    }else if (!identical(length(w_y), length(y))) {
      stop("The weight function does not return a numeric vector of the same length as y.")
    } else if (any(w_y < 0)) {
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
    } else if (any(diff(v_y) < -tol)) {
      message("The chaining function is not increasing.")
    }
  }
  
}
