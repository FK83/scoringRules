#' Default Weight and Chaining Functions
#' 
#' Get commonly used weight or chaining functions to use within weighted scoring rules.
#' The normal and logistic distribution, density, and survival functions are available.
#' Multivariate normal distribution functions are also available for multivariate scoring rules.
#' 
#' @param name name of the weight function to extract.
#' @param mu location parameter of the normal or logistic distribution.
#' @param sigma scale parameter of the normal or logistic distribution.
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
#' Allen, S., Ginsbourger, D. and J. Ziegel (2022): 
#' `Evaluating forecasts for high-impact events using transformed kernel scores', 
#' \emph{arXiv preprint} arXiv:2202.12732.
#' \doi{10.48550/arXiv.2202.12732}
#'  
#' @author Sam Allen
#' 
#' @seealso \link{scores_sample_univ_weighted} and \link{scores_sample_multiv_weighted} for weighted scoring rules.
#' 
#' @details 
#' Details will be added here
#' 
#' @examples
#' \dontrun{
#' Examples will be added here
#' }
#' 
#' @name get_weight_func
NULL

#' @rdname get_weight_func
#' @export
get_weight_func <- function (name = c("norm_cdf", "norm_surv", "norm_pdf", "logis_cdf", "logis_surv", "logis_pdf"),
                             mu = 0, sigma = 1, weight = TRUE) {
  input <- list(name = name, mu = mu, sigma = sigma, weight = weight)
  check_gwf_input(input)
  
  if (weight) {
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
  if (length(mu) > 1 || length(sigma) > 1) stop("mu and sigma must be single numeric values.")
  if (sigma <= 0) stop("sigma must be a single positive value.")
  
  if (!is.logical(weight)) stop("weight must be a logical.")
}
