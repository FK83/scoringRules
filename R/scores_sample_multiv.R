#' Multivariate Scoring Rules for Simulated Forecast Distributions
#' 
#' Compute multivariate scores of the form \eqn{S(y, dat)}, where \eqn{S} is a
#' proper scoring rule, \eqn{y} is a d-dimensional realization vector and 
#' \eqn{dat} is a simulated sample of multivariate forecasts. Available are the
#' energy score and the variogram score of order \eqn{p}.
#' 
#' @param y realized values (numeric vector of length d).
#' @param dat numeric matrix of data
#' (columns are simulation draws from multivariate forecast distribution).
#' @param w numeric matrix of weights for \code{dat} used in the variogram
#' score. If no weights are specified, constant weights with \eqn{w = 1}
#' are used.
#' @param p order of variogram score. Standard choices include \eqn{p = 1} and
#' \eqn{p = 0.5}.
#' 
#' @details
#' In the input matrix \code{dat} each column is expected to represent a sample
#' from the multivariate forecast distribution, the number of rows of \code{dat}
#' thus has to match the length of the observation vector \code{y}, and the
#' number of columns of \code{dat} is the number of simulated samples.
#' 
#' In \link{vs_sample} it is possible to specify a matrix \code{w} of
#' non-negative weights that allow to emphasize or downweight pairs of
#' component combinations based on subjective expert decisions. \code{w} is a
#' square matrix with dimensions equal to the length of the input vector
#' \code{y}, and the entry in the \eqn{i}-th row and \eqn{j}-th column of
#' \code{w} corresponds to the weight assigned to the combination of the
#' corresponding \eqn{i}-th and \eqn{j}-th component. For details and examples,
#' see Scheuerer and Hamill (2015).
#' 
#' @return
#' Value of the score. \emph{A lower score indicates a better forecast.}
#' 
#' @references
#' \emph{Energy score}
#' 
#' Gneiting, T., Stanberry, L.I., Grimit, E.P., Held, L. and
#' N.A. Johnson (2008):
#' `Assessing probabilistic forecasts of multivariate quantities, with an
#' application to ensemble predictions of surface winds',
#' TEST, 17, 211-235.
#' 
#' \emph{Variogram-based	proper scoring rules}
#' 
#' Scheuerer, M. and T.M. Hamill (2015):
#' `Variogram-based proper scoring rules for probabilistic forecasts of
#' multivariate quantities',
#' Monthly Weather Review, 143, 1321-1334.
#' 
#' @author Maximiliane Graeter, Sebastian Lerch, Fabian Krueger
#' 
#' @examples
#' d <- 10  # number of dimensions
#' m <- 50  # number of samples from multivariate forecast distribution
#' 
#' mu0 <- rep(0, d)
#' mu <- rep(1, d)
#' S0 <- S <- diag(d)
#' S[S==0] <- 0.1
#' S0[S0==0] <- 0.2
#' 
#' # generate samples from multivariate normal distributions
#' obs <- drop(mu0 + rnorm(d) \%*\% chol(S0))
#' fc_sample <- replicate(m, drop(mu + rnorm(d) \%*\% chol(S)))
#' 
#' es_sample(y = obs, dat = fc_sample)
#' 
#' # weighting matrix for variogram score
#' w_vs <- matrix(NA, nrow = d, ncol = d)
#' for(d1 in 1:d){for(d2 in 1:d){w_vs[d1,d2] <- 0.5^abs(d1-d2)}}
#' 
#' vs_sample(y = obs, dat = fc_sample) 
#' vs_sample(y = obs, dat = fc_sample, w = w_vs) 
#' vs_sample(y = obs, dat = fc_sample, w = w_vs, p = 1)
#' 
#' @name scores_sample_multiv
NULL

################################################################################
# energy score
#' @rdname scores_sample_multiv
#' @export
es_sample <- function(y, dat) {
  input <- list(y = y, dat = dat)
  check.multivsample(input)
  energyscoreC(y, dat)
}

################################################################################
# variogram score of order p
#' @rdname scores_sample_multiv
#' @export
vs_sample <- function(y, dat, w = NULL,  p = 0.5) {
  input <- list(y = y, dat = dat)
  check.multivsample(input)
  d <- length(y)
  
  # additional input checks for weighting matrix w and order p
  if (is.null(w)) {
    w <- matrix(1, nrow = d, ncol = d)
  } else {
    if (!is.matrix(w)) {
      stop("'w' is not a matrix ")
    }
    if (any(dim(w) != d)) {
      stop("Dimensions of 'w' do not fit")
    }
    if (any(w < 0)) {
      stop("Weighting matrix 'w' contains negative values")
    }
  }

  if (!is.numeric(p) || length(p) != 1 ){
    stop("Order 'p' must be numeric of length 1")
  } else if (p < 0) {
    stop("Order 'p' must be positive")
  }
  
  
  out <- 0
  for (i in 1:d) {
    for (j in 1:d) {
      vdat <- mean(abs(dat[i, ] - dat[j, ])^p)
      vy <- abs(y[i] - y[j])^p
      out <- out + w[i, j] * (vy - vdat)^2 
    }
  }
  
  return(out)
}

################################################################################
### input checks for multivariate scoring rules

check.multivsample <- function(input) {
  input_isnumeric <- sapply(input, is.numeric)
  if (!all(input_isnumeric)) {
    stop(paste("Non-numeric input:", 
               paste(names(input)[!input_isnumeric], collapse=", ")))
  }
  
  if (!is.vector(input$y)) {
    stop("'y' is not a vector")
  } 
  
  if (!is.matrix(input$dat)) {
    stop("'dat' is not a matrix ")
  }
  
  if (length(input$y) != dim(input$dat)[1]) {
    stop("Dimensions of 'y' and 'dat' do not fit")
  }
}
