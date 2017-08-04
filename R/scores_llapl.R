#' Calculating scores for the log-Laplace distribution
#'
#' @param y vector of observations.
#' @param locationlog vector of location parameters on the log scale.
#' @param scalelog vector of positive scale parameters on the log scale.
#' @return A vector of score values.
#' @name scores_llapl
NULL

#' @rdname scores_llapl
#' @export
crps_llapl <- function(y, locationlog, scalelog) {
  scalelog[scalelog < 0 | scalelog >= 1] <- NaN
  y1 <- y
  y1[y1 < 0] <- 0
  z <- (log(y1) - locationlog) / scalelog
  p <- 0.5 + 0.5 * sign(z) * pexp(abs(z))
  c1 <- y*(2*p - 1)
  c2 <- ifelse (z < 0,
                (1 - (2*p)^(1 + scalelog)) / (1 + scalelog),
                - (1 - (2*(1-p))^(1 - scalelog)) / (1 - scalelog)
  )
  c3 <- scalelog / (4 - scalelog^2)
  c1 + exp(locationlog) * (c2 + c3)
}

#' @rdname scores_llapl
#' @export
logs_llapl <- function(y, locationlog, scalelog) {
  -log(fllapl(y, locationlog, scalelog))
}
  


check_crps_llapl <- function(input) {
  required <- c("y", "locationlog", "scalelog")
  checkNames1(required, names(input))
  checkNumeric(input)
  checkVector(input)
  
  if (any(input$scalelog <= 0))
    stop("Parameter 'scalelog' contains non-positive values.")
  if (any(input$scalelog >= 1)) {
    stop(paste("Parameter 'scalelog' contains values greater or equal to 1.",
               "The CRPS does not exist."))
  }
}

check_logs_llapl <- function(input) {
  required <- c("y", "locationlog", "scalelog")
  checkNames1(required, names(input))
  checkNumeric(input)
  checkVector(input)
  
  if (any(input$scalelog <= 0))
    stop("Parameter 'scalelog' contains non-positive values.")
}
