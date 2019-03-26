#' Scoring Rules for a Vector of Moments
#' 
#' Calculate scores (DSS, ESS) given observations and moments of the predictive distributions.
#' 
#' @param y vector of realized values.
#' @param mean vector of mean values.
#' @param var vector of variance values.
#' @param skew vector of skewness values. 
#' 
#' @return
#' Value of the score. \emph{A lower score indicates a better forecast.}
#' 
#' @references
#' \emph{Dawid-Sebastiani score:}
#' 
#' Dawid, A.P. and P. Sebastiani (1999):
#' 'Coherent dispersion criteria for optimal experimental design'
#' The Annals of Statistics, 27, 65-81.
#'  
#' \emph{Error-spread score:}
#'  
#' Christensen, H.M., I.M. Moroz, and T.N. Palmer (2015):
#' `Evaluation of ensemble forecast uncertainty using a new proper score:
#' Application to medium-range and seasonal forecasts',
#' Quarterly Journal of the Royal Meteorological Society, 141, 538-549.
#'  
#' @author Alexander Jordan, Sebastian Lerch
#' 
#' @details 
#' The skewness of a random variable \eqn{X} is the third standardized moment 
#' \deqn{E[(\frac{X-\textnormal{mean}}{\sqrt{\textnormal{var}}})^3].} 
#' 
#' 
#' @name scores_moments
NULL

#' @rdname scores_moments
#' @export
dss_moments <- function(y, mean = 0, var = 1) {
  (y - mean)^2 / var + log(var)
}

#' @rdname scores_moments
#' @export
ess_moments <- function(y, mean = 0, var = 1, skew = 0) {
  e <- mean - y
  (var - e^2 - e*sqrt(var)*skew)^2
}
