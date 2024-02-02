#' Default Weight and Chaining Functions
#' 
#' Get commonly used weight or chaining functions to use within weighted scoring rules.
#' The normal and logistic distribution, density, and survival functions are available.
#' Multivariate normal distribution functions are also available for multivariate scoring rules.
#' 
#' @param name name of the weight function to extract.
#' @param mu location parameter(s) of the normal or logistic distribution.
#' @param sigma scale parameter(s) of the normal or logistic distribution.
#' @param weight logical specifying whether to return a weight function (\code{weight = TRUE})
#'  or chaining function (\code{weight = FALSE}).
#'  
#' @return
#' A weight or chaining function.
#' 
#' @references
#' 
#' Gneiting, T. and R. Ranjan (2011): 
#' `Comparing density forecasts using threshold-and quantile-weighted scoring rules', 
#' \emph{Journal of Business & Economic Statistics} 29, 411-422. 
#' \doi{10.1198/jbes.2010.08110}
#' 
#' Allen, S., Ginsbourger, D. and J. Ziegel (2023): 
#' `Evaluating forecasts for high-impact events using transformed kernel scores', 
#' \emph{SIAM/ASA Journal on Uncertainty Quantification} 11, 906-940.
#' \doi{10.1137/22M1532184}
#'  
#' @author Sam Allen
#' 
#' @seealso \link{scores_sample_univ_weighted} and \link{scores_sample_multiv_weighted} for weighted scoring rules.
#' 
#' @details 
#' The weighted scoring rules in \link{scores_sample_univ_weighted} and \link{scores_sample_multiv_weighted}
#' require a weight or chaining function argument (\code{weight_func} or \code{chain_func}) to target
#' particular outcomes. \code{get_weight_func()} can be used to obtain the relevant R function 
#' corresponding to commonly-used weight and chaining functions.
#'
#' These commonly-used weight and chaining functions correspond to cumulative distribution functions (cdf's),
#' probability density function (pdf's) and survival functions of the normal and logistic distributions.
#' The \code{name} argument specifies the desired weight or chaining function. This must be one of \code{'norm_cdf'},
#' \code{'norm_pdf'}, \code{'norm_surv'}, \code{'logis_cdf'}, \code{'logis_pdf'} and \code{'logis_surv'}, 
#' corresponding to the cdf, pdf and survival functions of the normal and logistic distribution, respectively.
#' 
#' \code{mu} and \code{sigma} represent the location and scale parameters of the normal or logistic distribution.
#' 
#' \code{weight} is a logical that specifies whether a weight or chaining function should be returned: 
#' if \code{weight = TRUE} (the default) the weight function is returned, and if \code{weight = FALSE}
#' the chaining function is returned.
#' 
#' The normal weight and chaining functions are applicable in both the univariate and multivariate setting.
#' In the univariate case, \code{mu} and \code{sigma} should be single numeric values. In the multivariate case, 
#' \code{'norm_cdf'} and \code{'norm_pdf'} represent the cdf and pdf of the multivariate normal distribution, 
#' with mean vector \code{mu} and covariance matrix \code{diag(sigma)}. Here, \code{mu} and \code{sigma} are
#' vectors with length equal to the dimension of the multivariate outcomes.
#' 
#' Note that \code{get_weight_func()} can currently only return multivariate weight and chaining 
#' functions corresponding to the multivariate normal distribution with a diagonal covariance matrix.
#' 
#' 
#' @examples
#' \dontrun{
#' 
#' ## univariate
#' # generate data
#' y <- rnorm(10)
#' sample_fc <- matrix(rnorm(100), nrow = 10)
#' 
#' # normal cdf
#' mu <- 1
#' sigma <- 1
#' 
#' weight_func <- get_weight_func("norm_cdf", mu = mu, sigma = sigma)
#' chain_func <- get_weight_func("norm_cdf", mu = mu, sigma = sigma, weight = FALSE)
#' owcrps_sample(y = y, dat = sample_fc, weight_func = weight_func)
#' twcrps_sample(y = y, dat = sample_fc, chain_func = chain_func)
#' 
#' # results are the same if the weight function is input manually
#' weight_func <- function(x) pnorm(x, mu, sigma)
#' chain_func <- function(x) (x - mu)*pnorm(x, mu, sigma) + (sigma^2)*dnorm(x, mu, sigma)
#' owcrps_sample(y = y, dat = sample_fc, weight_func = weight_func)
#' twcrps_sample(y = y, dat = sample_fc, chain_func = chain_func)
#' 
#' 
#' # logistic pdf
#' mu <- 0
#' sigma <- 1
#' 
#' weight_func <- get_weight_func("logis_pdf", mu = mu, sigma = sigma)
#' chain_func <- get_weight_func("logis_pdf", mu = mu, sigma = sigma, weight = FALSE)
#' owcrps_sample(y = y, dat = sample_fc, weight_func = weight_func)
#' twcrps_sample(y = y, dat = sample_fc, chain_func = chain_func)
#' 
#'
#' # normal survival function 
#' mu <- -1
#' sigma <- 1
#' 
#' weight_func <- get_weight_func("norm_surv", mu = mu, sigma = sigma)
#' chain_func <- get_weight_func("norm_surv", mu = mu, sigma = sigma, weight = FALSE)
#' owcrps_sample(y = y, dat = sample_fc, weight_func = weight_func)
#' twcrps_sample(y = y, dat = sample_fc, chain_func = chain_func)
#'
#'
#' ## multivariate
#' d <- 3  # number of dimensions
#' m <- 10  # number of samples from multivariate forecast distribution
#' 
#' # generate samples from multivariate normal distributions
#' mu0 <- rep(0, d)
#' mu <- rep(1, d)
#' S0 <- S <- diag(d)
#' S0[S0==0] <- 0.2
#' S[S==0] <- 0.1
#' 
#' y <- drop(mu0 + rnorm(d) %*% chol(S0))
#' sample_fc <- replicate(m, drop(mu + rnorm(d) %*% chol(S)))
#'
#' # component-wise normal cdf
#' mu <- rep(1, d)
#' sigma <- rep(1, d)
#' 
#' weight_func <- get_weight_func("norm_cdf", mu = mu, sigma = sigma)
#' chain_func <- get_weight_func("norm_cdf", mu = mu, sigma = sigma, weight = FALSE)
#' owes_sample(y = y, dat = sample_fc, weight_func = weight_func)
#' twes_sample(y = y, dat = sample_fc, chain_func = chain_func)
#'
#' }
#' 
#' @name get_weight_func
NULL

#' @rdname get_weight_func
#' @export
get_weight_func <- function (name = "norm_cdf", mu = 0, sigma = 1, weight = TRUE) {
  input <- list(name = name, mu = mu, sigma = sigma, weight = weight)
  check_gwf_input(input)
  
  if (weight) {
    if ((identical(length(mu), 1L) && identical(length(sigma), 1L))) {
      if (name == "norm_cdf") {
        weight_func <- function(x) pnorm(x, mean = mu, sd = sigma)
      } else if (name == "norm_surv") {
        weight_func <- function(x) 1 - pnorm(x, mean = mu, sd = sigma)
      } else if (name == "norm_pdf") {
        weight_func <- function(x) dnorm(x, mean = mu, sd = sigma)
      } else if (name == "logis_cdf") {
        weight_func <- function(x) plogis(x, location = mu, scale = sigma)
      } else if (name == "logis_surv") {
        weight_func <- function(x) 1 - plogis(x, location = mu, scale = sigma)
      } else if (name == "logis_pdf") {
        weight_func <- function(x) dlogis(x, location = mu, scale = sigma)
      }
    } else {
      if (name == "norm_cdf") {
        weight_func <- function(x) prod(pnorm(x, mean = mu, sd = sigma))
      } else if (name == "norm_surv") {
        weight_func <- function(x) prod(1 - pnorm(x, mean = mu, sd = sigma))
      } else if (name == "norm_pdf") {
        weight_func <- function(x) prod(dnorm(x, mean = mu, sd = sigma))
      }
    }
    return(weight_func)
  } else {
    if (name == "norm_cdf") {
      chain_func <- function(x) 
        (x - mu)*pnorm(x, mean = mu, sd = sigma) + (sigma^2)*dnorm(x, mean = mu, sd = sigma)
    } else if (name == "norm_surv") {
      chain_func <- function(x) 
        x - (x - mu)*pnorm(x, mean = mu, sd = sigma) - (sigma^2)*dnorm(x, mean = mu, sd = sigma)
    } else if (name == "norm_pdf") {
      chain_func <- function(x) pnorm(x, mean = mu, sd = sigma)
    } else if (name == "logis_cdf") {
      chain_func <- function(x) sigma*log(1 + exp((x - mu)/sigma))
    } else if (name == "logis_surv") {
      chain_func <- function(x) x - sigma*log(1 + exp((x - mu)/sigma))
    } else if (name == "logis_pdf") {
      chain_func <- function(x) plogis(x, location = mu, scale = sigma)
    }
    return(chain_func)
  }
}


################################################################################
# checks for the input to get_weight_func
check_gwf_input <- function(input) {
  name <- input$name
  mu <- input$mu
  sigma <- input$sigma
  weight <- input$weight
  
  if (!is.character(name)) {
    stop("name must be one of 'norm_cdf', 'norm_surv', 'norm_pdf', 'logis_cdf', 'logis_surv', 'logis_pdf'.")
  } else if (!(name %in% c('norm_cdf', 'norm_surv', 'norm_pdf', 'logis_cdf', 'logis_surv', 'logis_pdf'))) {
    stop("name must be one of 'norm_cdf', 'norm_surv', 'norm_pdf', 'logis_cdf', 'logis_surv', 'logis_pdf'.")
  }
  
  if (!is.numeric(mu)) stop("mu must be numeric.")
  if (!is.numeric(sigma)) stop("sigma must be numeric.")
  if (!identical(length(mu), length(sigma))) stop("mu and sigma must be the same length.")
  if (any(sigma <= 0)) stop("sigma must only contain positive values.")
  
  if (!is.logical(weight)) stop("weight must be a logical.")
  
  if (!(identical(length(mu), 1L) && identical(length(sigma), 1L)) && 
      name %in% c('logis_cdf', 'logis_surv', 'logis_pdf')) {
    stop("Logistic weight functions are not available for multivariate forecasts. Use instead one of 'norm_cdf', 'norm_surv' and 'norm_pdf'.")
  }
}
