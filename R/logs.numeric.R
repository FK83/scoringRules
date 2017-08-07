#' Logarithmic Score for Parametric Forecast Distributions
#' 
#' Calculate the logarithmic score (LogS) given observations
#' and parameters of a family of distributions.
#' 
#' @param y Vector of realized values.
#' @param family String which specifies the parametric family; current options:
#' \code{"2pexp", "2pnorm", "beta", "exp", "exp2",
#' "exponential", "gamma", "gev", "gpd", "lapl",
#' "laplace", "llapl", "llogis", "lnorm", "log-laplace", "log-logistic",
#' "log-normal", "logis", "logistic", "mixnorm", "mixture-normal", "nbinom",
#' "negative-binomial", "norm", "normal", "pois", "poisson", "t", "tlogis",
#' "tnorm", "tt", "two-piece-exponential", "two-piece-normal", "unif", "uniform"}.
#' @param ... Vectors of parameter values; expected input depends on the chosen
#' \code{family}. See details below.
#' 
#' @return Vector of score values.
#' \emph{A lower score indicates a better forecast.}
#' 
#' @author Alexander Jordan, Fabian Krueger, Sebastian Lerch
#' 
#' @details
#' The parameters supplied to each of the functions are numeric vectors:
#' \enumerate{
#'  \item Distributions defined on the real line:
#'    \itemize{
#'      \item
#'        \code{"laplace"} or \code{"lapl"}:
#'        \code{location} (real-valued location parameter),
#'        \code{scale} (positive scale parameter);
#'        see \code{\link{logs_lapl}}
#'      \item
#'        \code{"logistic"} or \code{"logis"}:
#'        \code{location} (real-valued location parameter),
#'        \code{scale} (positive scale parameter);
#'        see \code{\link{logs_logis}}
#'      \item
#'        \code{"normal"} or \code{"norm"}:
#'        \code{mean}, \code{sd} (mean and standard deviation);
#'        see \code{\link{logs_norm}}
#'      \item
#'        \code{"normal-mixture"} or \code{"mixture-normal"} or \code{"mixnorm"}:
#'        \code{m} (mean parameters),
#'        \code{s} (standard deviations),
#'        \code{w} (weights);
#'        see \code{\link{logs_mixnorm}};
#'        note: matrix-input for parameters
#'      \item
#'        \code{"t"}:
#'        \code{df} (degrees of freedom),
#'        \code{location} (real-valued location parameter),
#'        \code{scale} (positive scale parameter);
#'        see \code{\link{logs_t}}
#'      \item
#'        \code{"two-piece-exponential"} or \code{"2pexp"}:
#'        \code{location} (real-valued location parameter),
#'        \code{scale1}, \code{scale2} (positive scale parameters);
#'        see \code{\link{logs_2pexp}}
#'      \item
#'        \code{"two-piece-normal"} or \code{"2pnorm"}:
#'        \code{location} (real-valued location parameter),
#'        \code{scale1}, \code{scale2} (positive scale parameters);
#'        see \code{\link{logs_2pnorm}}
#'    }
#'  \item Distributions for non-negative random variables:
#'    \itemize{
#'      \item
#'        \code{"exponential"} or \code{"exp"}:
#'        \code{rate} (positive rate parameter);
#'        see \code{\link{logs_exp}}
#'      \item
#'        \code{"gamma"}:
#'        \code{shape} (positive shape parameter),
#'        \code{rate} (positive rate parameter),
#'        \code{scale} (alternative to \code{rate});
#'        see \code{\link{logs_gamma}}
#'      \item
#'        \code{"log-laplace"} or \code{"llapl"}:
#'        \code{locationlog} (real-valued location parameter),
#'        \code{scalelog} (positive scale parameter);
#'        see \code{\link{logs_llapl}}
#'      \item
#'        \code{"log-logistic"} or \code{"llogis"}:
#'        \code{locationlog} (real-valued location parameter),
#'        \code{scalelog} (positive scale parameter);
#'        see \code{\link{logs_llogis}}
#'      \item
#'        \code{"log-normal"} or \code{"lnorm"}:
#'        \code{locationlog} (real-valued location parameter),
#'        \code{scalelog} (positive scale parameter);
#'        see \code{\link{logs_lnorm}}
#'    }
#'  \item Distributions with flexible support and/or point masses:
#'    \itemize{
#'      \item
#'        \code{"beta"}:
#'        \code{shape1}, \code{shape2} (positive shape parameters),
#'        \code{lower}, \code{upper} (lower and upper limits);
#'        see \code{\link{logs_beta}}
#'      \item
#'        \code{"uniform"} or \code{"unif"}:
#'        \code{min}, \code{max} (lower and upper limits);
#'        see \code{\link{logs_unif}}
#'      \item
#'        \code{"exp2"}:
#'        \code{location} (real-valued location parameter),
#'        \code{scale} (positive scale parameter);
#'        see \code{\link{logs_exp2}}
#'      \item
#'        \code{"gev"}:
#'        \code{location} (real-valued location parameter),
#'        \code{scale} (positive scale parameter),
#'        \code{shape} (real-valued shape parameter);
#'        see \code{\link{logs_gev}}
#'      \item
#'        \code{"gpd"}:
#'        \code{location} (real-valued location parameter),
#'        \code{scale} (positive scale parameter),
#'        \code{shape} (real-valued shape parameter);
#'        see \code{\link{logs_gpd}}
#'      \item
#'        \code{"tlogis"}:
#'        \code{location} (location parameter),
#'        \code{scale} (scale parameter),
#'        \code{lower}, \code{upper} (lower and upper limits);
#'        see \code{\link{logs_tlogis}}
#'      \item
#'        \code{"tnorm"}:
#'        \code{location} (location parameter),
#'        \code{scale} (scale parameter),
#'        \code{lower}, \code{upper} (lower and upper limits);
#'        see \code{\link{logs_tnorm}}
#'      \item
#'        \code{"tt"}:
#'        \code{df} (degrees of freedom),
#'        \code{location} (location parameter),
#'        \code{scale} (scale parameter),
#'        \code{lower}, \code{upper} (lower and upper limits);
#'        see \code{\link{logs_tt}}
#'    }
#'  \item Distributions of discrete variables:
#'    \itemize{
#'      \item
#'        \code{"negative-binomial"} or \code{"nbinom"}:
#'        \code{size} (positive dispersion parameter),
#'        \code{prob} (success probability),
#'        \code{mu} (mean, alternative to \code{prob});
#'        see \code{\link{logs_nbinom}}
#'      \item
#'        \code{"poisson"} or \code{"pois"}:
#'        \code{lambda} (positive mean);
#'        see \code{\link{logs_pois}}
#'    }
#' }
#' All numerical arguments should be of the same length.
#' An exception are scalars of length 1, which will be recycled.
#' 
#' @examples
#' logs(y = 1, family = "normal", mean = 0, sd = 2) 
#' logs(y = rnorm(20), family = "normal", mean = 1:20, sd = sqrt(1:20))
#' 
#' ## Arguments can have different lengths:
#' logs(y = rnorm(20), family = "normal", mean = 0, sd = 2)
#' logs(y = 1, family = "normal", mean = 1:20, sd = sqrt(1:20))
#' 
#' ## Mixture of normal distributions requires matrix input for parameters:
#' mval <- matrix(rnorm(20*50), nrow = 20)
#' sdval <- matrix(runif(20*50, min = 0, max = 2), nrow = 20)
#' weights <- matrix(rep(1/50, 20*50), nrow = 20)
#' logs(y = rnorm(20), family = "mixnorm", m = mval, s = sdval, w = weights)
#' 
#' @seealso \code{\link{crps.numeric}}
#' 
#' @export logs.numeric
#' @export
logs.numeric <- function(y, family, ...) {
  family <- getFamily(family, "logs")
  checkInput <- get(paste0("check_logs_", family))
  calculateLogS <- get(paste0("logs_", family))
  
  input <- list(y = y, ...)
  checkInput(input)
  out <- do.call(calculateLogS, input)
  
  if (any(is.na(out))) {
    warning("Missing LogS values.")
  }
  
  out
}
