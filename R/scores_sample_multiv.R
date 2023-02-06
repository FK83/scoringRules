#' Multivariate Scoring Rules for Simulated Forecast Distributions
#' 
#' Compute multivariate scores of the form \eqn{S(y, dat)}, where \eqn{S} is a
#' proper scoring rule, \eqn{y} is a d-dimensional realization vector and 
#' \eqn{dat} is a simulated sample of multivariate forecasts. Three scores
#' are available: The energy score, a score based on a Gaussian kernel 
#' (\link{mmds_sample}, see details below) and the variogram score of order \eqn{p}.
#' 
#' @param y realized values (numeric vector of length d).
#' @param dat numeric matrix of data
#' (columns are simulation draws from multivariate forecast distribution).
#' @param w numeric vector of weights for forecast draws (length equal to number of columns of \code{dat})
#' @param w_vs numeric matrix of weights for \code{dat} used in the variogram
#' score. This matrix must be square and symmetric, with all elements being non-negative.
#' If no weights are specified, constant weights (with all elements of \code{w_vs} 
#' equal to one) are used.
#' @param p order of variogram score. Standard choices include \eqn{p = 1} and
#' \eqn{p = 0.5}.
#' 
#' @details
#' In the input matrix \code{dat} each column is expected to represent a sample
#' from the multivariate forecast distribution, the number of rows of \code{dat}
#' thus has to match the length of the observation vector \code{y}, and the
#' number of columns of \code{dat} is the number of simulated samples.
#' 
#' In \link{es_sample} and \link{mmds_sample} it is possible to specify a vector \code{w} of weights 
#' attached to each forecast draw (i.e. each column of matrix \code{dat}). These
#' weights are analogous to the input \code{w} in \link{crps_sample}. 
#' 
#' In \link{vs_sample} it is possible to specify a matrix \code{w_vs} of
#' non-negative weights that allow to emphasize or downweight pairs of
#' component combinations based on subjective expert decisions. \code{w_vs} is a
#' square, symmetric matrix with dimensions equal to the length of the input vector
#' \code{y}, and the entry in the \eqn{i}-th row and \eqn{j}-th column of
#' \code{w_vs} corresponds to the weight assigned to the combination of the
#' corresponding \eqn{i}-th and \eqn{j}-th component. A small example is provided below. 
#' For details and further examples, see Scheuerer and Hamill (2015).
#' 
#' The `MMD score' in \link{mmds_sample} is a kernel scoring rule as described in 
#' Gneiting and Raftery (2007, Section 5). As for all other scores, 
#' we use a negative orientation, such that a smaller score corresponds to a better
#' forecast. We use a Gaussian kernel with standard deviation 1. This kernel is
#' the same as the one obtained by setting \code{rbfdot(sigma = .5)} in the 
#' R package kernlab (Karatzoglou et al., 2004). The naming prefix `MMD' is 
#' motivated by the machine learning literature on two sample testing 
#' (e.g. Gretton et al., 2012), where this type of kernel function is popular. 
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
#' TEST, 17, 211-235. \doi{10.1007/s11749-008-0114-x}
#' 
#' \emph{MMD score}
#'  
#' Gneiting, T. and A.E. Raftery (2007):
#' `Strictly proper scoring rules, prediction and estimation',
#' Journal of the American Statistical Association 102, 359-378. \doi{10.1198/016214506000001437}
#' 
#' Gretton, A., Borgwardt, K. M., Rasch, M. J., Schölkopf, B. and
#' A. Smola (2012): `A kernel two-sample test', Journal of` Machine 
#' Learning Research, 13, 723-773.
#' 
#' Karatzoglou, A., Smola, A., Hornik, K. and Zeileis A. (2004). 
#' kernlab - An S4 Package for Kernel Methods in R. Journal of Statistical 
#' Software 11, 1-20. \doi{10.18637/jss.v011.i09} 
#' 
#' \emph{Variogram-based proper scoring rules}
#' 
#' Scheuerer, M. and T.M. Hamill (2015):
#' `Variogram-based proper scoring rules for probabilistic forecasts of
#' multivariate quantities',
#' Monthly Weather Review, 143, 1321-1334. \doi{10.1175/mwr-d-14-00269.1}
#' 
#' @author Maximiliane Graeter, Sebastian Lerch, Fabian Krueger
#' 
#' @seealso \code{\link{scores_sample_multiv_weighted}} for weighted versions of the scoring rules documented here.
#' 
#' @examples
#' d <- 10  # number of dimensions
#' m <- 50  # number of samples from multivariate forecast distribution
#' 
#' # parameters for multivariate normal example
#' mu0 <- rep(0, d)
#' mu <- rep(1, d)
#' S0 <- S <- diag(d)
#' S0[S0==0] <- 0.2
#' S[S==0] <- 0.1
#' 
#' # generate samples from multivariate normal distributions
#' obs <- drop(mu0 + rnorm(d) %*% chol(S0))
#' fc_sample <- replicate(m, drop(mu + rnorm(d) %*% chol(S)))
#' 
#' # compute Energy Score
#' es_sample(y = obs, dat = fc_sample)
#' 
#' # in the univariate case, Energy Score and CRPS are the same
#' # illustration: Evaluate forecast sample for the first variable
#' es_sample(y = obs[1], dat = fc_sample[1, , drop = FALSE])
#' crps_sample(y = obs[1], dat = fc_sample[1, ])
#' 
#' # illustration of observation weights for Energy Score
#' # example: equal weights for first half of draws; zero weights for other draws
#' w <- rep(c(1, 0), each = .5*m)/(.5*m)
#' es_sample(y = obs, dat = fc_sample, w = w)
#' 
#' # weighting matrix for variogram score
#' # note that, unlike for w, weights in w_vs refer to dimensions 
#' # (rows of dat) rather than draws (cols of dat)
#' w_vs <- outer(1:d, 1:d, function(x, y) .5^abs(x-y))
#' 
#' vs_sample(y = obs, dat = fc_sample) 
#' vs_sample(y = obs, dat = fc_sample, w_vs = w_vs) 
#' vs_sample(y = obs, dat = fc_sample, w_vs = w_vs, p = 1)
#' 
#' 
#' @name scores_sample_multiv
NULL

################################################################################
# energy score
#' @rdname scores_sample_multiv
#' @export
es_sample <- function(y, dat, w = NULL) {
  input <- list(y = y, dat = dat)
  check.multivsample(input)
  w <- w.helper.multiv(dat, w)
  es <- esC_xy(y, dat, w) - .5*esC_xx(dat, w)
  return(es)
}

################################################################################
# MMD score
#' @rdname scores_sample_multiv
#' @export
mmds_sample <- function(y, dat, w = NULL) {
  input <- list(y = y, dat = dat)
  check.multivsample(input)
  w <- w.helper.multiv(dat, w)
  # note that order of xx and xy parts is reverse to Energy Score
  # (since underlying kernels are in reverse orientation)
  mmds <- .5*mmdsC_xx(dat, w) - mmdsC_xy(y, dat, w)
  return(mmds)
}

################################################################################
# variogram score of order p
#' @rdname scores_sample_multiv
#' @export
vs_sample <- function(y, dat, w = NULL, w_vs = NULL, p = 0.5) {
  input <- list(y = y, dat = dat)
  check.multivsample(input)
  d <- length(y)
  # additional input checks for weighting matrix w_vs and order p
  if (!is.null(w_vs)) {
    if (!is.matrix(w_vs)) {
      stop("'w_vs' is not a matrix ")
    }
    if (any(dim(w_vs) != d)) {
      stop("Dimensions of 'w_vs' do not fit")
    }
    if (any(w_vs < 0)) {
      stop("Weighting matrix 'w_vs' contains negative values")
    }
    if (!isSymmetric(w_vs)) {
      stop("Weighting matrix 'w_vs' is not symmetric")
    }
  }
  if (!is.numeric(p) || length(p) != 1 ){
    stop("Order 'p' must be numeric of length 1")
  } else if (p < 0) {
    stop("Order 'p' must be positive")
  }
  # Compute score 
  if (is.null(w)) {
    if (is.null(w_vs)) {
      out <- vsC(y, dat, p)
    } else {
      out <- vsC_w_vs(y, dat, w_vs, p)
    }
  } else {
    w_vs <- matrix(1, nrow = length(y), ncol = length(y))
    w <- w.helper.multiv(dat, w)
    out <- vsC_w(y, dat, w_vs, w, p)
  }
  return(out)
}

################################################################################
# helper functions
# input checks for multivariate scoring rules
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
# function to check or generate weights in `es_sample' and `mmds_sample` 
w.helper.multiv <- function(dat, w){
  m <- ncol(dat)
  if (is.null(w)){
    w_final <- rep(1, m)/m
  } else {
    if (any(w < 0) || (length(w) != m) || !is.vector(w)){
      stop("w is of wrong format")
    } else {
      sum_w <- sum(w)
      if (abs(sum_w - 1) > 1e-5){
        message("Weights in w don't add to one - they have been re-scaled to one")
      }
      w_final <- w/sum_w
    }
  }
  return(w_final)
}