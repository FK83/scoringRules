#' Scoring Rules for Simulated Forecast Distributions
#' 
#' Compute scores of the form \eqn{S(y, dat)}, where \eqn{S} is a proper
#' scoring rule, \eqn{y} is a scalar realization and \eqn{dat} is a simulated
#' sample of data.
#' 
#' @param y realized value (numeric vector of length one).
#' @param dat Numeric vector of data
#'  (simulation draws from forecast distribution).
#' @param method String, approximation method used in \link{crps_sample}.
#'  Currently implemented: "edf" (empirical distribution function) and
#'  "kde" (kernel density estimation).
#' @param w Numeric vector of weights for \code{dat} (if \code{method = "edf"}).
#' @param bw Numeric scalar for bandwidth (if \code{method = "kde"}).
#' @param num_int Logical; if TRUE numerical integration is used
#'  (if \code{method = "kde"}).
#' @param show_messages Logical; on/off switch for messages (does not affect
#'  warnings and error, which are always printed).
#'  
#' @return
#' Value of the score. \emph{A lower score indicates a better forecast.}
#' 
#' @references
#' \emph{Evaluating simulation based forecast distributions:}
#' 
#' Krueger, F., S. Lerch, T.L. Thorarinsdottir and T. Gneiting (2017),
#' `Probabilistic forecasting and comparative model assessment based on
#' Markov Chain Monte Carlo output', working paper,
#' Heidelberg Institute for Theoretical Studies,
#' available at \url{http://arxiv.org/abs/1608.06802}.
#'  
#' \emph{Empirical quantile decomposition of the CRPS:}
#'  
#' Laio, F. and S. Tamea (2007),
#' `Verification tools for probabilistic forecasts of continuous
#' hydrological variables',
#' Hydrology and Earth System Sciences, 11, 1267-1277.
#'  
#' @author Alexander Jordan, Fabian Krueger
#'  
#' @details
#' When using \code{method = "edf"}, \link{crps_sample} employs an empirical
#' version of the quantile decomposition of the CRPS (Laio and Tamea, 2007).
#' When using \code{method = "kde"}, the function does kernel density
#' estimation using a Gaussian kernel. The bandwidth (\code{bw}) can be
#' specified manually, in which case it must be a positive number. If
#' \code{bandwidth == NULL}, the bandwidth is selected using the core function
#' \link{bw.nrd}. Numerical integration may speed up computation for
#' \link{crps_sample} in case of large samples \code{dat} if
#' \code{method = "kde"}.
#'  
#' @examples
#' sample <- rnorm(1e4)
#' obs <- 1
#' crps_sample(y = obs, dat = sample)  
#' logs_sample(y = obs, dat = sample)
#' 
#' ## exact solutions using analytical closed form expression:
#' crps(y = obs, family = "normal", mean = 0, sd = 1)
#' logs(y = obs, family = "normal", mean = 0, sd = 1)
#' 
#' @name scores_sample_univ
#' @importFrom stats bw.nrd
NULL

#' @rdname scores_sample_univ
#' @export
crps_sample <- function(y, dat, method = "edf", w = NULL, bw = NULL, 
                        num_int = FALSE, show_messages = TRUE) {
  input <- list(y = y, dat = dat)
  if (!is.null(w)) {
    input$w <- w
  }
  if (!is.null(bw)) {
    input$bw <- bw
  }
  check.sample(input)
  
  # Further checks for inconsistent/not implemented inputs
  if (method == "edf"){
    if (!is.null(bw)){
      warning("Parameter 'bw' is ignored for edf method.")
    }
    if (num_int){
      warning("Parameter 'num_int' is ignored for edf method.")
    }
  }
  if (method == "kde") {
    if (!is.null(w)) {
      warning("Parameter 'w' is ignored for kde method.")
    }
  }
  if (!method %in% c("edf", "kde")){
    stop("Unexpected choice of method - please select either edf or kde")
  }
  
  # Score computations
  if (method == "edf") {
    out <- crps.edf(dat = dat, y = y, w = w)
  } else if (method == "kde") {
    # Bandwidth
    if (is.null(bw)) {
      bw <- bw.nrd(dat)
    }
    if (num_int == FALSE) {
      # Exact formula
      out <- crps.kdens(dat = dat, y = y, bw = bw)
    } else {
      # Numerical integration
      FF <- function(x) {
        sapply(x, function(zz) mean(pnorm(zz, mean = dat, sd = bw)))
      }
      aux1 <- integrate(function(s) FF(s)^2,-Inf, y, rel.tol = 1e-6)$value
      aux2 <- integrate(function(s) (1 - FF(s))^2, y, Inf, rel.tol = 1e-6)$value
      out <- aux1 + aux2
      # Message
      if (show_messages) 
        message("Used numerical integration for CRPS computation (tolerance = 1e-6).")
    }
  } 
  return(out)
}


#' @rdname scores_sample_univ
#' @export
logs_sample <- function(y, dat, bw = NULL, show_messages = TRUE) {
  input <- list(y = y, dat = dat)
  if (!is.null(bw)) {
    input$bw <- bw
  }
  check.sample(input)
  
  if (show_messages)
    message("Using the log score with kernel density estimation tends to be fragile -- see KLTG (2016) for details.")
  if (is.null(bw)) bw <- bw.nrd(dat)
  w <- rep(1 / length(dat), length(dat))
  s <- rep(bw, length(dat))
  out <- lsmixnC(w = w, m = dat, s = s, y = y)
  return(out)
}
