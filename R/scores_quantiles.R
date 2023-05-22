#' Quantile and interval scores
#' @param y vector of observations
#' @param x vector of quantile predictions
#' @param x_lower,x_upper vector of quantile predictions (lower and upper endpoints of prediction intervals)
#' @param alpha quantile level of interest 
#' @param target_coverage target (i.e., nominal) coverage level of prediction interval
#' @return A vector of score values. Smaller values indicate better forecasts. Note that 
#' the interval score refers to the central prediction interval at level \code{target_coverage}.
#' @references 
#' Quantile score:
#' Koenker, R. and G. Bassett (1978): `Regression quantiles', Econometrica 46, 33-50. \doi{https://doi.org/10.2307/1913643}
#' Interval score:
#' Gneiting, T. and A.E. Raftery (2007):
#' `Strictly proper scoring rules, prediction and estimation',
#' Journal of the American Statistical Association 102, 359-378. \doi{10.1198/016214506000001437}
#' @examples
#' qs(y = 1, x = qnorm(.1), alpha = .1)
#' qs(y = 1, x = qnorm(.9), alpha = .9)
#' # Interval score is proportional to sum of two quantile scores
#' ints(y = 1, x_lower = qnorm(.1), x_upper = qnorm(.9), target_coverage = .8)
#' @name scores_quantiles
#' @export
qs <- function(y, x, alpha){
  # input checks
  input <- list(y = y, x = x)
  checkNumeric(input)
  checkVector(input)
  if (!is.numeric(alpha) || length(alpha) != 1 ){
    stop("alpha must be numeric of length 1")
  } else if (alpha <= 0 || alpha >= 1) {
    stop("alpha must satisfy 0 < alpha < 1")
  }
  # score calculation
  ((y < x)-alpha)*(x-y)
}

#' interval score
#' @export
#' @rdname scores_quantiles
ints <- function(y, x_lower, x_upper, target_coverage){
  # input checks
  input <- list(y = y, x_lower = x_lower, x_upper = x_upper)
  checkNumeric(input)
  checkVector(input)
  if (!is.numeric(target_coverage) || length(target_coverage) != 1 ){
    stop("target_coverage must be numeric of length 1")
  } else if (target_coverage <= 0 || target_coverage >= 1) {
    stop("target_coverage must satisfy 0 < target_coverage < 1")
  }
  if (any(input$x_lower > input$x_upper)){
    stop("'x_lower' contains values greater than corresponding values in 'x_upper'.")  
  }
  # get interval score as sum of two quantile scores
  aux1 <- qs(y, x_lower, .5*(1-target_coverage)) 
  aux2 <- qs(y, x_upper, .5*(1+target_coverage)) 
  (2/(1-target_coverage))*(aux1 + aux2)
}