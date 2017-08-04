#' Calculating scores for the log-normal distribution
#'
#' @param y vector of observations.
#' @param meanlog an alternative way to specify \code{locationlog}.
#' @param sdlog an alternative way to specify \code{scalelog}.
#' @param locationlog vector of location parameters on the log scale.
#' @param scalelog vector of positive scale parameters on the log scale.
#' @return A vector of score values.
#' @name scores_lnorm
#' @importFrom stats plnorm dlnorm
NULL

#' @rdname scores_lnorm
#' @export
crps_lnorm <- function(y, meanlog = 0, sdlog = 1,
                       locationlog = meanlog, scalelog = sdlog) {
  if (!missing(meanlog) && !missing(locationlog))
    stop("specify 'meanlog' or 'locationlog' but not both")
  if (!missing(sdlog) && !missing(scalelog))
    stop("specify 'sdlog' or 'scalelog' but not both")
  c1 <- y * (2 * plnorm(y, locationlog, scalelog) - 1)
  c2 <- 2 * exp(locationlog + 0.5 * scalelog^2)
  c3 <- plnorm(y, locationlog + scalelog^2, scalelog) +
    pnorm(scalelog / sqrt(2)) - 1
  c1 - c2 * c3
}


#' @rdname scores_lnorm
#' @export
logs_lnorm <- function(y, meanlog = 0, sdlog = 1,
                       locationlog = meanlog, scalelog = sdlog) {
  if (!missing(meanlog) && !missing(locationlog))
    stop("specify 'meanlog' or 'locationlog' but not both")
  if (!missing(sdlog) && !missing(scalelog))
    stop("specify 'sdlog' or 'scalelog' but not both")
  -dlnorm(y, locationlog, scalelog, log = TRUE)
}
  

check_crps_lnorm <- function(input) {
  required <- list(c("y", "meanlog", "sdlog"),
                   c("y", "locationlog", "scalelog"))
  checkNames2(required, names(input))
  checkNumeric(input)
  checkVector(input)
  
  if ("sdlog" %in% names(input)) {
    if (any(input$sdlog <= 0))
      stop("Parameter 'sdlog' contains non-positive values.")
  }
  if ("scalelog" %in% names(input)) {
    if (any(input$scalelog <= 0))
      stop("Parameter 'scalelog' contains non-positive values.")
  }
}

check_logs_lnorm <- check_crps_lnorm
