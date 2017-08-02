#' Generic Scoring Rule Calculation
#' 
#' @description
#' Generic functions for calculating the Continuous Ranked Probability Score
#' and the Logarithmic Score of R objects.
#' 
#' \code{scoringRules} provides default methods
#' (\link{crps.numeric}, \link{logs.numeric}) to calculate scores of forecasts
#' that are members of families of parametric distributions.
#' 
#' @param y an object for which the score is to be calculated
#' @param ... further arguments passed to or from other methods
#' 
#' @return Returns a vector of scores. One for each forecast-observation pair.
#' 
#' @seealso
#' \link{crps.numeric}, \link{logs.numeric}
#' 
#' @name scores
NULL


#' @rdname scores
#' @export
crps <- function(y, ...) UseMethod("crps")

#' @rdname scores
#' @export
logs <- function(y, ...) UseMethod("logs")
