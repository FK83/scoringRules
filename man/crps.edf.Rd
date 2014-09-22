\name{crps}
\alias{crps.edf}
\alias{crps.int}
\alias{crps.krep}
\alias{crps.kdens}
\alias{crps.mixn}
\alias{crps.t}
\title{
Functions to compute the crps for a forecast distribution and realization
}
\description{
The crps (continuous ranked probability score) is a performance measure for distribution forecasts, defined for a cumulative distribution function \eqn{F} and realizing observation \eqn{y}. The functions documented here compute the crps for various forms of \eqn{F}.

\link{crps.edf}, \link{crps.krep} and \link{crps.kdens} take a sample of data drawn from \eqn{F} as an input (this form of \eqn{F} is common in MCMC). The three functions differ with respect to the estimation method used (see details below). \link{crps.int} is suitable when \eqn{F} is available as a function. \link{crps.mixn} and \link{crps.t} are for parametric forms of \eqn{F} (mixture of normals and t, respectively).
}
\usage{
crps.edf(dat, y, w = NULL)
crps.int(F, y)
crps.krep(dat, y)
crps.mixn(m, v, y, exact = TRUE, w = NULL)
crps.kdens(dat, y, exact = TRUE, bw = NULL)
crps.t(m, v, df, y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
Numeric vector of data (simulation draws from forecast distribution).
}
  \item{w}{\link{crps.edf}: Numeric vector of weights for \code{dat}. \link{crps.mixn}: Numeric vector of weights for mixture components. In both cases: Elements must be between zero and one and sum to one; default is equal weights.}
\item{F}{R function \code{F(z)}, returning the cdf of the forecast distribution at a vector \code{z}}
\item{m, v}{\link{crps.mixn}: Vectors of means and variances for all mixture components. Lengths must be equal, and equal to the length of \code{w}. \link{crps.t}: Mean and variance parameters, each of length one.}
\item{exact}{logical, set to TRUE for exact analytical formula, and to FALSE for much quicker approximation using numerical integration. Applies to \link{crps.mixn} and \link{crps.kdens}.}
\item{df}{Degrees of freedom parameter for the t distribution (numeric vector of length one).}
\item{y}{Realizing observation (numeric vector of length one).}
\item{bw}{Bandwidth parameter supplied to \link{crps.kdens} (numeric vector of length one, positive). Defaults to using \link{bw.SJ} for choosing the parameter.}
}

\value{
Value of the CRPS. 
}
\references{
General background on scoring rules and probabilistic forecasting: 

Gneiting, T. and A. Raftery (2007): ```Strictly Proper Scoring Rules, Prediction and Estimation'', Journal of the American Statistical Association 102, 359-378.

Gneiting, T. and M. Katzfuss (2014): ``Probabilistic Forecasting'', Annual Review of Statistics and Its Application 1, 125-151. 

Evaluating simulation based forecast distributions:

Krueger, F., S. Lerch, T. Thorarinsdottir and T. Gneiting (2015), ``Evaluating MCMC-based Predictive Distributions'', working paper, Heidelberg Institute for Theoretical Studies, \emph{in preparation}.

Overview on analytical and numerical ways to evaluate the crps:

Jordan, A. (2015): ``Evaluating the CRPS'', working paper, Heidelberg Institute for Theoretical Studies, \emph{in preparation}.

Paper which first proposed the crps:

Matheson, J.E. and R.L. Winkler (1976), ``Scoring Rules for Continuous Probability Distributions'', {Management Science} 22, 1087-1096.

}
\author{
Fabian Krueger
}
\details{
\link{crps.krep} is the empirical analogue of the kernel representation (Eq. 21 in Gneiting and Raftery, 2007); however, it performs a permutation of the draws in \eqn{dat} to correct for possible autocorrelation. \link{crps.edf} is based on the algorithm by Hersbach (2000, "Decomposition of the Continuous Ranked Probability Score for Ensemble Prediction Systems", Weather and Forecasting). \link{crps.kdens} is based on kernel density estimation, using a Gaussian kernel. If \code{exact == TRUE}, the latter function is quite slow when \code{dat} is large. The quick approximation obtained when setting \code{exact == FALSE} is usually very accurate.
}

\examples{
\dontrun{

# Example where all functions should yield very similar output

set.seed(287)
n <- 10000
y <- 1
dat <- rnorm(n)

crps.edf(dat, y)
crps.int(pnorm, y)
crps.krep(dat, y)
crps.mixn(0, 1, y)
crps.t(0, 1, 5000, y)
crps.kdens(dat, y)

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Scoring Rules}
\keyword{crps}
