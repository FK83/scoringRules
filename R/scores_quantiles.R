#' Quantile and interval scores
#' @description Compute quantile and interval scores, for  
#' given quantile predictions 
#' @param y vector of observations
#' @param x vector of quantile predictions
#' @param dat vector or matrix (depending on \code{y}; see details)
#'  of simulation draws from forecast distribution. 
#' @param x_lower,x_upper vector of quantile predictions (lower and upper endpoints of prediction intervals)
#' @param alpha quantile level of interest 
#' @param target_coverage target (i.e., nominal) coverage level of prediction interval
#' @param type integer between 1 and 9 that is passed on to stats function \link[stats]{quantile} (specifies algorithm for 
#' empirical quantile estimation; defaults to 7)
#' @param show_messages logical; display of messages (does not affect warnings and errors).
#' @param w vector of observation weights (currently ignored)
#' @return A vector of score values. Smaller values indicate better forecasts. Note that 
#' the interval score refers to the central prediction interval at level \code{target_coverage}.
#' @references 
#' Quantile score
#' 
#' Koenker, R. and G. Bassett (1978): `Regression quantiles', Econometrica 46, 33-50. \doi{https://doi.org/10.2307/1913643}
#' 
#' Interval score
#' 
#' Gneiting, T. and A.E. Raftery (2007):
#' `Strictly proper scoring rules, prediction and estimation',
#' Journal of the American Statistical Association 102, 359-378. \doi{10.1198/016214506000001437}
#' 
#' @seealso The syntax of \code{\link{qs_sample}} and \code{\link{ints_sample}} is analogous to the functions documented on \code{\link{scores_sample_univ}}.
#'  
#' @details 
#' For a vector \code{y} of length n, \code{dat} should be given as a matrix
#' with n rows. If \code{y} has length 1, then \code{dat} may be a vector.
#' @name scores_quantiles
#' @examples
#' # Example 1: Illustrate that interval score is proportional to sum of two quantile scores
#' target_coverage <- .8
#' # corresponding quantile levels
#' alpha_1 <- .5*(1-target_coverage)
#' alpha_2 <- 1-.5*(1-target_coverage)
#' # compute interval score
#' ints_quantiles(y = 1, x_lower = qnorm(alpha_1), 
#' x_upper = qnorm(alpha_2), target_coverage = target_coverage)
#' # compute sum of quantile scores (scaled by 2/(1-target_coverage))
#' (2/(1-target_coverage))*(qs_quantiles(y = 1, x = qnorm(alpha_1), alpha = alpha_1) + 
#' qs_quantiles(y = 1, x = qnorm(alpha_2), alpha = alpha_2))
#' 
#' # Example 2: Compare exact to simulated quantile forecast from standard normal distribution
#' qs_quantiles(y = 1, x = qnorm(.1), alpha = .1)
#' qs_sample(y = 1, dat = rnorm(500), alpha = .1)
#' @export
qs_quantiles <- function(y, x, alpha){
  # input checks
  input <- list(y = y, x = x)
  checkNumeric(input)
  checkVector(input)
  check_alpha(alpha)
  # score calculation
  ((y < x)-alpha)*(x-y)
}

# function to check provided quantile level
check_alpha <- function(alpha){
  if (!is.numeric(alpha) || !identical(length(alpha), 1L)){
    stop("alpha must be numeric of length 1")
  } else if (isTRUE(alpha <= 0) || isTRUE(alpha >= 1)) {
    stop("alpha must satisfy 0 < alpha < 1")
  }
}

# function to check inputs for interval score
check_ints <- function(target_coverage, 
                       x_lower, x_upper){
  if (!is.numeric(target_coverage) || !identical(length(target_coverage), 1L)){
    stop("target_coverage must be numeric of length 1")
  } else if (isTRUE(target_coverage <= 0) || isTRUE(target_coverage >= 1)) {
    stop("target_coverage must satisfy 0 < target_coverage < 1")
  }
  if (isTRUE(any(x_lower > x_upper))){
    stop("'x_lower' contains values greater than corresponding values in 'x_upper'.")  
  }  
}

# interval score
#' @export
#' @rdname scores_quantiles
ints_quantiles <- function(y, x_lower, x_upper, target_coverage){
  # input checks
  input <- list(y = y, x_lower = x_lower, x_upper = x_upper)
  checkNumeric(input)
  checkVector(input)
  check_ints(target_coverage, x_lower, x_upper)
  # get interval score as sum of two quantile scores
  aux1 <- qs_quantiles(y, x_lower, .5*(1-target_coverage)) 
  aux2 <- qs_quantiles(y, x_upper, .5*(1+target_coverage)) 
  (2/(1-target_coverage))*(aux1 + aux2)
}

qs_edf <- function(y, dat, alpha, w = NULL, type = 7) {
  # compute empirical alpha quantile 
  q_alpha <- quantile(dat, probs = alpha, names = FALSE, type = type)
  f <- function(s){
    qs_quantiles(y = s, x = q_alpha, alpha = alpha)
  }
  sapply(y, f)
}

#' @rdname scores_quantiles
#' @export
qs_sample <- function(y, dat, alpha, w = NULL, 
                      type = 7, show_messages = TRUE) {
  input <- list(y = y, dat = dat)
  if ( (!is.null(w)) && isTRUE(show_messages) ){
    message("Parameter 'w' is currently ignored for qs_sample.")
    w <- NULL 
  }
  if (identical(length(y), 1L) && is.vector(dat)) {
    check_sample(input)
    qs_edf(y = y, dat = dat, alpha = alpha, type = type)
  } else {
    check_sample2(input)
    sapply(seq_along(y),
           function(i) qs_edf(y = y[i], dat = dat[i, ], alpha = alpha, 
                              type = type))
  }
}

ints_edf <- function(y, dat, target_coverage, w = NULL, type = 7) {
  # compute relevant empirical quantiles
  l <- c(.5*(1-target_coverage), 1-.5*(1-target_coverage))
  q_alpha <- quantile(dat, probs = l, names = FALSE, type = type)
  f <- function(s){
    ints_quantiles(y = s, x_lower = q_alpha[1], x_upper = q_alpha[2], 
                   target_coverage = target_coverage)
  }
  sapply(y, f)
}

#' @rdname scores_quantiles
#' @export
ints_sample <- function(y, dat, target_coverage, w = NULL, 
                        type = 7, show_messages = TRUE) {
  input <- list(y = y, dat = dat)
  if ( (!is.null(w)) && isTRUE(show_messages) ){
    message("Parameter 'w' is currently ignored for ints_sample.")
    w <- NULL 
  }
  if (identical(length(y), 1L) && is.vector(dat)) {
    check_sample(input)
    ints_edf(y = y, dat = dat, target_coverage = target_coverage, 
             type = type)
  } else {
    check_sample2(input)
    sapply(seq_along(y),
           function(i) ints_edf(y = y[i], dat = dat[i, ], 
                                target_coverage = target_coverage, 
                                type = type))
  }
}