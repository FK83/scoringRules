\name{Log Score}
\alias{ls.mixn}
\alias{ls.kdens}
\title{
Functions to compute the logarithmic score for a forecast distribution and realization
}
\description{
The logarithmic (or ``log'') score is a performance measure for distribution forecasts, defined for a probability density function \eqn{f} and realizing observation \eqn{y}. The functions documented here compute the log score for various forms of \eqn{f}.
}
\usage{
ls.mixn(m, v, y, w = NULL)
ls.kdens(dat, y, bw = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
Numeric vector of data (simulation draws from forecast distribution).
}
  \item{w}{\link{crps.edf}: Numeric vector of weights for \code{dat}. \link{crps.mixn}: Numeric vector of weights for mixture components. In both cases: Elements must be between zero and one and sum to one; default is equal weights.}
\item{m, v}{Vectors of means and variances for all mixture components. Lengths must be equal, and equal to the length of \code{w}.}
\item{y}{Realizing observation (numeric vector of length one).}
\item{bw}{Bandwidth parameter supplied to \link{ls.kdens} (numeric vector of length one, positive). Defaults to using \link{bw.SJ} for choosing the parameter.}
}

\value{
Value of the log score. 
}
\references{
General background on scoring rules and probabilistic forecasting: 

Gneiting, T. and A. Raftery (2007): ```Strictly Proper Scoring Rules, Prediction and Estimation'', Journal of the American Statistical Association 102, 359-378.

Gneiting, T. and M. Katzfuss (2014): ``Probabilistic Forecasting'', Annual Review of Statistics and Its Application 1, 125-151. 

Evaluating simulation based forecast distributions:

Krueger, F., S. Lerch, T. Thorarinsdottir and T. Gneiting (2015), ``Evaluating MCMC-based Predictive Distributions'', working paper, Heidelberg Institute for Theoretical Studies, \emph{in preparation}.

Paper which first proposed the log score:

Good, I.J. (1952): ``Rational Decisions'', Journal of the Royal Statistical Society (B) 14, 107-114.

}
\author{
Fabian Krueger
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Scoring Rules}
\keyword{log score}
