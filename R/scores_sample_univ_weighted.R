#' Weighted Scoring Rules for Simulated Forecast Distributions
#' 
#' Calculate weighted scores (threshold-weighted CRPS, outcome-weighted CRPS, censored likelihood score) given observations and draws from the predictive distributions.
#' 
#' @param y vector of realized values.
#' @param dat vector or matrix (depending on \code{y}; see details)
#'  of simulation draws from forecast distribution. 
#' @param a optional; numeric lower bound for indicator weight function \code{w(z) = 1{a < z < b}}.
#' @param b optional; numeric upper bound for indicator weight function \code{w(z) = 1{a < z < b}}.
#' @param method string; approximation method. Options:
#'  "edf" (empirical distribution function) and
#'  "kde" (kernel density estimation).
#' @param w optional; vector or matrix (matching \code{dat}) of weights for method \code{"edf"}.
#' @param bw optional; vector (matching \code{y}) of bandwidths for kernel density
#' estimation; see details.
#' @param num_int logical; if TRUE numerical integration is used for method \code{"kde"}.
#' @param show_messages logical; display of messages (does not affect
#'  warnings and errors).
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
#' \emph{Outcome-weighted CRPS:}
#'  
#' Holzmann, H. and B. Klar (2017):
#' `Focusing on regions of interest in forecast evaluation',
#' \emph{Annals of Applied Statistics} 11, 2404-2431. 
#' \doi{10.1214/17-AOAS1088}
#' 
#' \emph{Censored likelihood score:}
#' 
#' Diks, C., Panchenko, V. and D. Van Dijk (2011):
#' `Likelihood-based scoring rules for comparing density forecasts in tails',
#' \emph{Journal of Econometrics} 163, 215-230.
#' \doi{10.1016/j.jeconom.2011.04.001}
#'  
#' @author Sam Allen
#'  
#' @details 
#' For a vector \code{y} of length n, \code{dat} should be given as a matrix
#' with n rows. If \code{y} has length 1, then \code{dat} may be a vector.
#' 
#' \code{\link{twcrps_sample}} transforms \code{y} and \code{dat} according to the
#' weight function, and then employs \code{\link{crps_sample}}. See the documentation
#' for \code{\link{crps_sample}} for further details.
#' 
#' The \code{w} argument is used to weight the draws from the predictive distribution. 
#' This argument is also present in the unweighted scores (e.g. \code{\link{crps_sample}}).
#' This does not weight particular outcomes within the weighted scoring rules.  
#' 
#' @name scores_sample_univ_weighted
NULL

#' @rdname scores_sample_univ_weighted
#' @export
twcrps_sample <- function(y, dat, a = -Inf, b = Inf, method = "edf", w = NULL, bw = NULL, 
                          num_int = FALSE, show_messages = TRUE) {
  input <- list(lower = a, upper = b)
  check_weight(input)
  if (method == "edf") {
    v_y <- pmin(pmax(y, a), b)
    v_dat <- pmin(pmax(dat, a), b)
    crps_sample(y = v_y, dat = v_dat, method, w, bw, num_int, show_messages)
  } else if (method == "kde") {
    stop("twcrps_sample is not yet available for method kde")
  } else {
    stop("Unexpected choice of method - please select either 'edf' or 'kde'.")
  }
  
}


#' @rdname scores_sample_univ_weighted
#' @export
owcrps_sample <- function(y, dat, a = -Inf, b = Inf, method = "edf", w = NULL, bw = NULL, 
                          num_int = FALSE, show_messages = TRUE, comp = TRUE) {
  input <- list(lower = a, upper = b)
  check_weight(input)
  weight_func <- function(x) as.numeric(x > a & x < b)
  w_y <- sapply(y, weight_func)
  if (method == "edf"){
    if (identical(length(y), 1L) && is.vector(dat)) {
      w_dat <- sapply(dat, weight_func)
    } else {
      w_dat <- array(sapply(dat, weight_func), dim(dat))
    }
    if (is.null(w)) {
      w <- w_dat
    }else {
      w <- w*w_dat
    }
    score <- crps_sample(y, dat, method, w = w, bw, num_int, show_messages)*w_y
    # complement the owCRPS with the Brier score
    if (comp) {
      if (identical(length(y), 1L) && is.vector(dat)) {
        int_wf <- mean(dat < b) - mean(dat < a)
      }else {
        int_wf <- rowMeans(dat < b) - rowMeans(dat < a)
      }
      score <- score + w_y*(1 - int_wf)^2 + (1 - w_y)*(int_wf^2)
    }
  }
  if (method == "kde") {
    stop("owcrps_sample is not yet available for method kde")
  } else {
    stop("Unexpected choice of method - please select either 'edf' or 'kde'.")
  }
  return (score)
}


#' @rdname scores_sample_univ_weighted
#' @export
clogs_sample <- function(y, dat, a = -Inf, b = Inf, w = NULL, bw = NULL, 
                          show_messages = TRUE, comp = TRUE) {
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
  score <- w_y*(ls + log(int_wf))
  if (comp) {
    score <- score - w_y*log(int_wf) - (1 - w_y)*log(1 - int_wf)
  }
  return (score)
}


#### input checks for weighted scoring rules ####
check_weight <- function(input) {
  a <- input$lower
  b <- input$upper
  
  if (length(a) > 1 | length(b) > 1) {
    stop("The lower and upper bounds in the weight function must be single numeric values.")
  } else if (a >= b) {
    stop("The lower bound in the weight function is not smaller than the upper bound.")
  }
  
}
