#' Generic Scoring Rule Calculation
#' 
#' @description
#' Generic functions for calculating the Continuous Ranked Probability Score
#' and the Logarithmic Score of R objects.
#' 
#' \code{scoringRules} provides default methods
#' (\code{\link{crps.numeric}}, \code{\link{logs.numeric}}) to calculate scores of forecasts
#' that are members of families of parametric distributions.
#' 
#' @param y an object for which the score is to be calculated
#' @param ... further arguments passed to or from other methods
#' 
#' @return Returns a vector of scores. One for each forecast-observation pair.
#' 
#' @details
#' The mean logarithmic score corresponds to the negative of the
#' log-likelihood \code{\link{logLik}}.
#' 
#' @references
#' \emph{General background and further references on scoring rules:}
#' 
#' Gneiting, T. and A.E. Raftery (2007):
#' `Strictly proper scoring rules, prediction and estimation',
#' Journal of the American Statistical Association 102, 359-378. \doi{10.1198/016214506000001437}
#' 
#' Gneiting, T. and M. Katzfuss (2014):
#' `Probabilistic forecasting',
#' Annual Review of Statistics and Its Application 1, 125-151. \doi{10.1146/annurev-statistics-062713-085831}
#' 
#' @seealso
#' \code{\link{crps.numeric}}, \code{\link{logs.numeric}}
#' 
#' @name scores
NULL


#' @rdname scores
#' @export
crps <- function(y, ...) UseMethod("crps")

#' @rdname scores
#' @export
logs <- function(y, ...) UseMethod("logs")
