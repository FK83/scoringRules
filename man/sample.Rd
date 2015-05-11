\name{sample}
\alias{crps.sample}
\alias{dss.sample}
\alias{ls.sample}
\alias{qs.sample}
\title{
Scoring Rules for Simulated Forecast Distributions
}
\description{
Compute scores of the form \eqn{S(dat, y)}, where \eqn{S} is a proper scoring rule, \eqn{dat} is a simulated sample of data, and \eqn{y} is a scalar realization. 
}
\usage{
crps.sample(dat, y, method = "edf", w = NULL, bw = NULL)
dss.sample(dat, y)
ls.sample(dat, y, bw = NULL)
qs.sample(dat, y, bw = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{dat}{Numeric vector of data (simulation draws from forecast distribution).}
\item{y}{Realized value (numeric vector of length one)}
\item{method}{String, approximation method used in \link{crps.sample}. Currently implemented: "edf" (empirical distribution function) and "kde" (kernel density estimation).}
\item{w}{Numeric vector of observation weights for \code{dat}.}
\item{bw}{Bandwidth parameter for kernel density estimator (numeric vector of length one)}
}

\value{Value of the score. \emph{A lower score indicates a better forecast.}}

\references{

\emph{Evaluating simulation based forecast distributions:}

Krueger, F., S. Lerch, T. Thorarinsdottir and T. Gneiting (2015), ``Evaluating MCMC-based Predictive Distributions'', working paper, Heidelberg Institute for Theoretical Studies, \emph{in preparation}.

}
\author{
Alexander Jordan, Fabian Krueger
}
\details{
When using \code{method = "edf"}, \link{crps.sample} employs the efficient algorithm by Hersbach (2000, "Decomposition of the Continuous Ranked Probability Score for Ensemble Prediction Systems", Weather and Forecasting). When using \code{meth = "kde"}, the function does kernel density estimation using a Gaussian kernel. The bandwidth (\code{bw}) can be specified manually, in which case it must be a positive number. If \code{bandwidth == NULL}, the bandwidth is selected using the core function \link{bw.SJ}.

\link{ls.sample} and \link{qs.sample} \emph{always} use kernel density estimation with a Gaussian kernel (the EDF option is not available since the two scores require a density). The choice of bandwidth is handled as in \link{crps.sample}.

\link{dss.sample} simply uses the empirical mean and variances of the observations in \code{dat}. 

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{crps}
\keyword{dss}
\keyword{qs}
\keyword{ls}
