#' Logarithmic Score for Parametric Forecast Distributions
#' 
#' Calculate the logarithmic score (LogS) given observations
#' and parameters of a family of distributions.
#' 
#' @param y Vector of realized values.
#' @param family String which specifies the parametric family; currently
#' implemented: "beta", "exponential", "gamma", "gev", "gpd", "laplace",
#' "log-laplace", "log-logistic", "log-normal", "logistic", "mixture-normal",
#' "negative-binomial", "normal", "poisson", "t", "two-piece-normal", "uniform".
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
#'        see \link{flapl}
#'      \item
#'        \code{"logistic"} or \code{"logis"}:
#'        \code{location} (real-valued location parameter),
#'        \code{scale} (positive scale parameter);
#'        see \link{Logistic}
#'      \item
#'        \code{"normal"} or \code{"norm"}:
#'        \code{mean}, \code{sd} (mean and standard deviation);
#'        see \link{Normal}
#'      \item
#'        \code{"t"}:
#'        \code{location} (real-valued location parameter),
#'        \code{scale} (positive scale parameter),
#'        \code{df} (degrees of freedom);
#'        see \link{ft}
#'      \item
#'        \code{"normal-mixture"} or \code{"mixnorm"}:
#'        \code{m} (mean parameters),
#'        \code{s} (standard deviations),
#'        \code{w} (weights);
#'        see \link{fmixnorm};
#'        note: matrix-input for parameters
#'      \item
#'        \code{"two-piece-exponential"} or \code{"2pexp"}:
#'        \code{location} (real-valued location parameter),
#'        \code{scale1}, \code{scale2} (positive scale parameters);
#'        see \link{f2pexp}
#'      \item
#'        \code{"two-piece-normal"} or \code{"2pnorm"}:
#'        \code{location} (real-valued location parameter),
#'        \code{scale1}, \code{scale2} (positive scale parameters);
#'        see \link{f2pnorm}
#'    }
#'  \item Distributions for non-negative random variables:
#'    \itemize{
#'      \item
#'        \code{"exponential"} or \code{"exp"}:
#'        \code{rate} (positive rate parameter);
#'        see \link{Exponential}
#'      \item
#'        \code{"gamma"}:
#'        \code{shape} (positive shape parameter),
#'        \code{rate} (positive rate parameter),
#'        \code{scale} (alternative to \code{rate});
#'        see \link{GammaDist}
#'      \item
#'        \code{"log-laplace"} or \code{"llapl"}:
#'        \code{locationlog} (real-valued location parameter),
#'        \code{scalelog} (positive scale parameter);
#'        see \link{fllapl}
#'      \item
#'        \code{"log-logistic"} or \code{"llogis"}:
#'        \code{locationlog} (real-valued location parameter),
#'        \code{scalelog} (positive scale parameter);
#'        see \link{fllogis}
#'      \item
#'        \code{"log-normal"} or \code{"lnorm"}:
#'        \code{locationlog} (real-valued location parameter),
#'        \code{scalelog} (positive scale parameter);
#'        see \link{Lognormal}
#'    }
#'  \item Distributions  for random variables with variable support:
#'    \itemize{
#'      \item
#'        \code{"normal"} or \code{"norm"}:
#'        \code{location} (location parameter),
#'        \code{scale} (scale parameter),
#'        \code{lower} (real-valued truncation parameter, lower bound),
#'        \code{upper} (real-valued truncation parameter, upper bound),
#'        \code{lmass} (point mass in lower bound, string "cens" or "trunc"),
#'        \code{umass} (point mass in upper bound, string "cens" or "trunc");
#'        see \link{fnorm}
#'      \item
#'        \code{"t"}:
#'        \code{location} (location parameter),
#'        \code{scale} (scale parameter),
#'        \code{df} (degrees of freedom),
#'        \code{lower} (real-valued truncation parameter, lower bound),
#'        \code{upper} (real-valued truncation parameter, upper bound),
#'        \code{lmass} (point mass in lower bound, string "cens" or "trunc"),
#'        \code{umass} (point mass in upper bound, string "cens" or "trunc");
#'        see \link{ft}
#'      \item
#'        \code{"logistic"} or \code{"logis"}:
#'        \code{location} (location parameter),
#'        \code{scale} (scale parameter),
#'        \code{lower} (real-valued truncation parameter, lower bound),
#'        \code{upper} (real-valued truncation parameter, upper bound),
#'        \code{lmass} (point mass in lower bound, string "cens" or "trunc"),
#'        \code{umass} (point mass in upper bound, string "cens" or "trunc");
#'        see \link{flogis} 
#'      \item
#'        \code{"exponential"} or \code{"exp"}:
#'        \code{location} (real-valued location parameter),
#'        \code{scale} (positive scale parameter),
#'        \code{mass} (point mass in \code{location});
#'        see \link{fexp}
#'      \item
#'        \code{"gpd"}:
#'        \code{location} (real-valued location parameter),
#'        \code{scale} (positive scale parameter),
#'        \code{shape} (real-valued shape parameter),
#'        \code{mass} (point mass in \code{location});
#'        see \link{fgpd}
#'      \item
#'        \code{"gev"}:
#'        \code{location} (real-valued location parameter),
#'        \code{scale} (positive scale parameter),
#'        \code{shape} (real-valued shape parameter);
#'        see \link{fgev}
#'    }
#'  \item Distribution for random variables defined on bounded intervals:
#'    \itemize{
#'      \item
#'        \code{"uniform"} or \code{"unif"}:
#'        \code{min}, \code{max} (lower and upper boundaries),
#'        \code{lmass}, \code{umass} (point mass in lower or upper boundary);
#'        see \link{Uniform}
#'      \item
#'        \code{"beta"}:
#'        \code{shape1}, \code{shape2} (positive parameters);
#'        see \link{Beta}
#'    }
#'  \item Distributions for random variables with discrete / infinite support:
#'    \itemize{
#'      \item
#'        \code{"poisson"} or \code{"pois"}:
#'        \code{lambda} (positive mean);
#'        see \link{Poisson}
#'      \item
#'        \code{"negative-binomial"} or \code{"nbinom"}:
#'        \code{size} (positive dispersion parameter),
#'        \code{prob} (success probability),
#'        \code{mu} (mean, alternative to \code{prob});
#'        see \link{NegBinomial}
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
#' @seealso \link{crps.numeric}
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
