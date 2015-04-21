\name{parametric}
\alias{crps.parametric}
\alias{dss.parametric}
\alias{ls.parametric}
\alias{qs.parametric}
\title{
Scoring Rules for Parametric Forecast Distributions
}
\description{
The functions documented here compute scores of the form \eqn{S(y, F)}, where \eqn{S} is a proper scoring rule, \eqn{y} is a vector of realizations, and \eqn{F} belongs to a parametric family of distributions. 
}
\usage{
crps.parametric(y, family, ...)
dss.parametric(y, family, ...)
ls.parametric(y, family, ...)
qs.parametric(y, family, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{y}{Vector of realized values}
\item{family}{String which specifies the parametric family; currently implemented: "beta", "exponential", "gamma", "gev", "gpd", "log-laplace", "log-logistic", "log-normal", "logistic", "mixture-normal", "negative-binomial", "normal", "poisson", "t", "truncated-normal", "two-piece-normal", "uniform".}
\item{...}{Vectors of parameter values; expected input depends on chosen \code{family}. See details below.}

}

\value{Vector of values of the score. \emph{A lower score indicates a better forecast.}}
\references{

\emph{General background and further references on scoring rules:}

Gneiting, T. and A. Raftery (2007): ```Strictly Proper Scoring Rules, Prediction and Estimation'', Journal of the American Statistical Association 102, 359-378.

Gneiting, T. and M. Katzfuss (2014): ``Probabilistic Forecasting'', Annual Review of Statistics and Its Application 1, 125-151. 

}
\author{
Alexander Jordan, Fabian Krueger, Sebastian Lerch
}
\details{
The parameters supplied to each of the functions are numeric vectors:
\enumerate{
\item Distributions for random variables with discrete / infinite support:
    \itemize{
      \item If \code{family == "poisson"}: \code{lambda} denotes the vector of positive means
      \item If \code{family == "negative-binomial"}: 
    }
\item Distribution for random variables defined on bounded intervals:
    \itemize{
      \item If \code{family == "beta"}: \code{shape1}, \code{shape2} ... simple reference to R help ? ...
      \item If \code{family == "uniform"}: \code{min} and \code{max} denote lower and upper limits of the distribution.
    }
\item Distributions defined on the real line:
    \itemize{
      \item If \code{family == "laplace"}: \code{location} and \code{scale} denote the location and scale parameter of the Laplace distribution with CDF ...
      \item If \code{family == "logistic"}: \code{location} and \code{scale} denote the location and scale parameter of the logistic distribution with CDF ...
      \item If \code{family == "normal"}: \code{mean} and \code{sd} denote the mean and standard deviation of the distribution.
      \item If \code{family == "normal-mixture"}: Matrices \code{m}, \code{s}, and \code{w}, where the \eqn{i}-th column of each matrix corresponds to the vector of mean values (\code{m[,i]}), standard deviations (\code{s[,i]}) and weights (\code{w[,i]}) of mixture component \eqn{i}. The weights in every row of \code{w} must sum to 1.
      \item If \code{family == "two-piece-normal"}: \code{m}, \code{s1} and \code{s2} denote the parameters of the two-piece normal distribution. For details, see Box A of Wallis (2004, "An Assessment of Bank of England and National Institute Inflation Forecast Uncertainties", National Institute Economic Review).
      \item If \code{family == "t"}:  \code{location}, \code{scale} (spread parameter) and \code{df} (degrees of freedom). Note that the parametrization is such that the implied variance of a t distributed variable is \eqn{scale^2 \times \frac{df}{df-2}}. 
    }
\item Distributions for non-negative random variables:
    \itemize{
      \item If \code{family == "exponential"}: \code{rate} denotes a vector of rate parameters
      \item If \code{family == "gamma"}:
      \item If \code{family == "log-laplace"}:
      \item If \code{family == "log-logistic"}:
      \item If \code{family == "log-normal"}: \code{meanlog} and \code{sdlog} denote location and scale parameter of the log-normal distribution.
    }
\item Distributions  for random variables with variable support:
    \itemize{
      \item If \code{family == "truncated-normal"}: \code{lb} (truncation parameter, lower bound), \code{m} (location parameter) and \code{s} (scale parameter). \code{m} and \code{s} correspond to mean and standard deviation of the normal distribution before the left-sided truncation in \code{lb}, respectively. [Reference to Gneiting and Thorarinsdottir (2010)?]
      \item If \code{family == "gpd"}: \code{loc} (location), \code{sc} (scale) and \code{sh} (shape) denote the parameters of the generalized Pareto distribution. For details, see Friederichs, F. and T.L. Thorarinsdottir (2012, "Forecast verification for extreme value distributions with an application to probabilistic peak wind prediction", Environmetrics 23, 579-594).
      \item If \code{family == "gev"}: \code{loc} (location), \code{sc} (scale) and \code{sh} (shape) denote the parameters of the generalized extreme value distribution. For details, see Friederichs, F. and T.L. Thorarinsdottir (2012, "Forecast verification for extreme value distributions with an application to probabilistic peak wind prediction", Environmetrics 23, 579-594).   
    }
}
}

\examples{
# Forecast distribution is N(0, 1), three different realizations (1, 2, and 3)
crps.parametric(c(1, 2, 3), family = "normal", mean = 0, sd = 1)

# Forecast distribution is a mixture of N(0, 1) and N(0, 3), single realization (2)
m <- matrix(0, ncol = 2)
s <- matrix(c(1, 3), ncol = 2)
w <- matrix(0.5, ncol = 2)
crps.parametric(2, family = "normal-mixture", m = m, s = s, w = w)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{crps}
\keyword{dss}
\keyword{qs}
\keyword{ls}
