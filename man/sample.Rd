\name{sample}
\alias{crps_sample}
\alias{logs_sample}
\title{
Scoring Rules for Simulated Forecast Distributions
}
\description{
Compute scores of the form \eqn{S(y, dat)}, where \eqn{S} is a proper scoring rule, \eqn{y} is a scalar realization and \eqn{dat} is a simulated sample of data. 
}
\usage{
crps_sample(y, dat, method = "edf", w = NULL, bw = NULL, 
            num_int = FALSE, show_messages = TRUE)
logs_sample(y, dat, bw = NULL, show_messages = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{y}{Realized value (numeric vector of length one).}
\item{dat}{Numeric vector of data (simulation draws from forecast distribution).}
\item{method}{String, approximation method used in \link{crps_sample}. Currently implemented: "edf" (empirical distribution function) and "kde" (kernel density estimation).}
\item{w}{Numeric vector of observation weights for \code{dat}.}
\item{bw}{Bandwidth parameter for kernel density estimator (numeric vector of length one)}
\item{num_int}{Logical of length one, whether or not to use numerical integration (may speed up computation for \link{crps_sample} and \link{qs_sample} in case of large samples \code{dat}).}
\item{show_messages}{Logical of length one, on/off switch for messages (does not affect warnings and errors, which are always printed)}
}

\value{Value of the score. \emph{A lower score indicates a better forecast.}}

\references{

\emph{Evaluating simulation based forecast distributions:}

Krueger, F., S. Lerch, T.L. Thorarinsdottir and T. Gneiting (2016), ``Evaluating MCMC-based Predictive Distributions'', working paper, Heidelberg Institute for Theoretical Studies, \emph{in preparation}.

}
\author{
Alexander Jordan, Fabian Krueger
}
\details{
When using \code{method = "edf"}, \link{crps_sample} employs an empirical version of the quantile decomposition of the CRPS (Laio, F. and Tamea, S. (2007), "Verification tools for probabilistic forecasts of continuous hydrological variables", Hydrology and Earth System Sciences, 11,
1267–1277). When using \code{method = "kde"}, the function does kernel density estimation using a Gaussian kernel. The bandwidth (\code{bw}) can be specified manually, in which case it must be a positive number. If \code{bandwidth == NULL}, the bandwidth is selected using the core function \link{bw.nrd}.

\link{logs_sample} \emph{always} uses kernel density estimation with a Gaussian kernel (the EDF option is not available since the score requires a density). The choice of bandwidth is handled as in \link{crps_sample}.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{crps}
\keyword{qs}
\keyword{logs}
