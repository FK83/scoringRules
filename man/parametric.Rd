\name{parametric}
\alias{crps.parametric}
\alias{dss.parametric}
\alias{ls.parametric}
\alias{qs.parametric}
\title{
Scoring Rules for Parametric Forecast Distributions
}
\description{
The functions documented here compute scores of the form \eqn{S(F, y)}, where \eqn{S} is a proper scoring rule, \eqn{F} belongs to a parametric family of distributions, and \eqn{y} is a scalar realization. 
}
\usage{
crps.parametric(family, parameters, y)
dss.parametric(family, parameters, y)
ls.parametric(family, parameters, y)
qs.parametric(family, parameters, y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{family}{String which specifies the parametric family; currently implemented: "normal", "truncated-normal", "two-piece-normal", "mixture-normal", "log-normal", "t", "gev", "gpd".}
\item{parameters}{List of parameters; expected format depends on chosen \code{family}. See details below.}
\item{y}{Realized value (numeric vector of length one)}
}

\value{Value of the score. \emph{A lower score indicates a better forecast.}}
\references{

\emph{General background and further references on scoring rules:}

Gneiting, T. and A. Raftery (2007): ```Strictly Proper Scoring Rules, Prediction and Estimation'', Journal of the American Statistical Association 102, 359-378.

Gneiting, T. and M. Katzfuss (2014): ``Probabilistic Forecasting'', Annual Review of Statistics and Its Application 1, 125-151. 

}
\author{
Alexander Jordan, Fabian Krueger, Sebastian Lerch
}
\details{
The parameters supplied to each of the functions are a list whose (named) elements are numeric vectors:
\enumerate{
\item If \code{family == "normal"}: List elements \code{m} and \code{s} denote the mean and standard deviation (std) of the distribution, respectively. Both are vectors of length one. 
\item If \code{family == "truncated-normal"}: List elements \code{lb} (truncation parameter, lower bound), \code{m} (location parameter) and \code{s} (scale parameter). All three are vectors of length one. \code{m} and \code{s} correspond to mean and standard deviation of the normal distribution before the left-sided truncation in \code{lb}, respectively. [Reference to Gneiting and Thorarinsdottir (2010) or similar?]
\item If \code{family == "two-piece-normal"}: List elements \code{m}, \code{s1} and \code{s2} (all of length one) denote the parameters of the two-piece normal distribution. For details, see Box A of Wallis (2004, "An Assessment of Bank of England and National Institute Inflation Forecast Uncertainties", National Institute Economic Review). 
\item If \code{family == "mixture-normal"}: Elements must be named \code{m}, \code{s}, and \code{w}. The \eqn{i}-th entry of each vector contains the mean, standard deviation and weight of mixture component \eqn{i}. All three vectors are of the same length (equal to the number of mixture components). 
\item If \code{family == "log-normal"}: List elements \code{loc} (location parameter) and \code{sh} (shape parameter). Both vectors of length one. 
\item If \code{family == "t"}: Elements \code{m} (mean), \code{s} (spread parameter) and \code{df} (degrees of freedom). All three vectors are of length one. Note that the parametrization is such that the implied variance of a t distributed variable is \eqn{s^2 \times \frac{df}{df-2}}. 
\item If \code{family == "gev"}: List elements \code{loc} (location), \code{sc} (scale) and \code{sh} (shape) denote the parameters of the generalized extreme value distribution. All three vectors are of length one. For details, see Friederichs, F. and T.L. Thorarinsdottir (2012, "Forecast verification for extreme value distributions with an application to probabilistic peak wind prediction", Environmetrics 23, 579-594).
\item If \code{family == "gpd"}: List elements \code{loc} (location), \code{sc} (scale) and \code{sh} (shape) denote the parameters of the generalized Pareto distribution. All three vectors are of length one. For details, see Friederichs, F. and T.L. Thorarinsdottir (2012, "Forecast verification for extreme value distributions with an application to probabilistic peak wind prediction", Environmetrics 23, 579-594).
}

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{crps}
\keyword{dss}
\keyword{qs}
\keyword{ls}
