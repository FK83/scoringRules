\name{parametric}
\alias{crps}
\alias{logs}
\keyword{logs}
\keyword{crps}
\title{
Scoring Rules for Parametric Forecast Distributions
}
\description{
Compute scores of the form \eqn{S(y, F)}, where \eqn{S} is a proper scoring rule, \eqn{y} is a vector of realizations, and \eqn{F} belongs to a parametric family of distributions. 
}
\usage{
crps(y, family, ...)
logs(y, family, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{y}{Vector of realized values}
\item{family}{String which specifies the parametric family; currently implemented: "beta", "exponential", "gamma", "gev", "gpd", "laplace", "log-laplace", "log-logistic", "log-normal", "logistic", "mixture-normal", "negative-binomial", "normal", "poisson", "t", "two-piece-normal", "uniform".}
\item{...}{Vectors of parameter values; expected input depends on chosen \code{family}. See details below.}

}

\value{Vector of values of the score. \emph{A lower score indicates a better forecast.}}
\references{

\emph{General background and further references on scoring rules:}

Gneiting, T. and A. Raftery (2007): `Strictly proper scoring rules, prediction and estimation', Journal of the American Statistical Association 102, 359-378.

Gneiting, T. and M. Katzfuss (2014): `Probabilistic forecasting', Annual Review of Statistics and Its Application 1, 125-151. 

\emph{Closed form expressions of the CRPS for specific distributions:}

Baran, S. and S. Lerch (2015): `Log-normal distribution based Ensemble Model Output Statistics models for probabilistic wind-speed forecasting', Quarterly Journal of the Royal Meteorological Society 141, 2289-2299. \emph{(Lognormal)}

Friederichs, P. and T.L. Thorarinsdottir (2012): `Forecast verification for extreme value distributions with an application to probabilistic peak wind prediction', Environmetrics 23, 579-594. \emph{(Generalized Extreme Value, Generalized Pareto)}

Gneiting, T., Larson, K., Westvelt III, A.H. and T. Goldman (2005): `Calibrated probabilistic forecasting using ensemble model output statistics and minimum CRPS estimation', Monthly Weather Review 133, 1098-1118. \emph{(Normal)}

Gneiting, T., Larson, K., Westrick, K., Genton, M.G. and E. Aldrich (2006): `Calibrated probabilistic forecasting at the stateline wind energy center: The regime-switching space-time method', Journal of the American Statistical Association 101, 968-979. \emph{(Censored normal)}

Gneiting, T. and T.L. Thorarinsdottir (2010): `Predicting inflation: Professional experts versus no-change forecasts', arXiv preprint arXiv:1010.2318. \emph{(Two-piece normal)}

Grimit, E.P.,  Gneiting, T., Berrocal, V.J. and N.A. Johnson (2006): `The continuous ranked probability score for circular variables and its application to mesoscale forecast ensemble verification', Quarterly Journal of the Royal Meteorological Society 132, 2925-2942. \emph{(Mixture of normals)}

Scheuerer, M. and D. Moeller (2015): `Probabilistic wind speed forecasting on a grid based on ensemble model output statistics', Annals of Applied Statistics 9, 1328-1349. \emph{(Gamma)}

Thorarinsdottir, T.L. and T. Gneiting (2010): `Probabilistic forecasts of wind speed: ensemble model output statistics by using heteroscedastic censored regression', Journal of the Royal Statistical Society (Series A) 173, 371-388. \emph{(Truncated normal)}

Wei, W. and L. Held (2014): `Calibration tests for count data', TEST 23, 787-205. \emph{(Poisson, Negative Binomial)}

\emph{Independent listing of closed-form solutions for the CRPS:}

Taillardat, M., Mestre, O., Zamo, M. and P. Naveau (2016): `Calibrated ensemble forecasts using quantile regression forests and ensemble model output statistics', Monthly Weather Review, in press. 

}
\author{
Alexander Jordan, Fabian Krueger, Sebastian Lerch
}
\details{
The parameters supplied to each of the functions are numeric vectors:
\enumerate{
\item Distributions defined on the real line:
    \itemize{
      \item \code{"laplace"} or \code{"lapl"}: \code{location} (real-valued location parameter), \code{scale} (positive scale parameter); see \link{flapl}
      \item \code{"logistic"} or \code{"logis"}: \code{location} (real-valued location parameter), \code{scale} (positive scale parameter); see \link{Logistic}
      \item \code{"normal"} or \code{"norm"}: \code{mean}, \code{sd} (mean and standard deviation); see \link{Normal}
      \item \code{"t"}:  \code{location} (real-valued location parameter), \code{scale} (positive scale parameter), \code{df} (degrees of freedom); see \link{ft}
      }
      \itemize{
      \item \code{"normal-mixture"} or \code{"mixnorm"}: \code{m} (mean parameters), \code{s} (standard deviations), \code{w} (weights); see \link{fmixnorm}; note: matrix-input for parameters
      \item \code{"two-piece-exponential"} or \code{"2pexp"}: \code{location} (real-valued location parameter), \code{scale1}, \code{scale2} (positive scale parameters); see \link{f2pexp}
      \item \code{"two-piece-normal"} or \code{"2pnorm"}: \code{location} (real-valued location parameter), \code{scale1}, \code{scale2} (positive scale parameters); see \link{f2pnorm}
    }
\item Distributions for non-negative random variables:
    \itemize{
      \item \code{"exponential"} or \code{"exp"}: \code{rate} (positive rate parameter); see \link{Exponential}
      \item \code{"gamma"}: \code{shape} (positive shape parameter), \code{rate} (positive rate parameter), \code{scale} (alternative to \code{rate}); see \link{GammaDist}
      \item \code{"log-laplace"} or \code{"llapl"}: \code{locationlog} (real-valued location parameter), \code{scalelog} (positive scale parameter); see \link{fllapl}
      \item \code{"log-logistic"} or \code{"llogis"}: \code{locationlog} (real-valued location parameter), \code{scalelog} (positive scale parameter); see \link{fllogis}
      \item \code{"log-normal"} or \code{"lnorm"}: \code{locationlog} (real-valued location parameter), \code{scalelog} (positive scale parameter); see \link{Lognormal}
    }
\item Distributions  for random variables with variable support:
    \itemize{
      \item \code{"normal"} or \code{"norm"}: \code{location} (location parameter), \code{scale} (scale parameter), \code{lower} (real-valued truncation parameter, lower bound), \code{upper} (real-valued truncation parameter, upper bound), \code{lmass} (point mass in lower bound, string "cens" or "trunc"), \code{umass} (point mass in upper bound, string "cens" or "trunc"); see \link{fnorm}
      \item \code{"t"}: \code{location} (location parameter), \code{scale} (scale parameter), \code{df} (degrees of freedom), \code{lower} (real-valued truncation parameter, lower bound), \code{upper} (real-valued truncation parameter, upper bound), \code{lmass} (point mass in lower bound, string "cens" or "trunc"), \code{umass} (point mass in upper bound, string "cens" or "trunc"); see \link{ft}
      \item \code{"logistic"} or \code{"logis"}: \code{location} (location parameter), \code{scale} (scale parameter), \code{lower} (real-valued truncation parameter, lower bound), \code{upper} (real-valued truncation parameter, upper bound), \code{lmass} (point mass in lower bound, string "cens" or "trunc"), \code{umass} (point mass in upper bound, string "cens" or "trunc"); see \link{flogis} 
      \item \code{"exponential"} or \code{"exp"}: \code{location} (real-valued location parameter), \code{scale} (positive scale parameter), \code{mass} (point mass in \code{location}); see \link{fexp}
      \item \code{"gpd"}: \code{location} (real-valued location parameter), \code{scale} (positive scale parameter), \code{shape} (real-valued shape parameter), \code{mass} (point mass in \code{location}); see \link{fgpd}
      \item \code{"gev"}: \code{location} (real-valued location parameter), \code{scale} (positive scale parameter), \code{shape} (real-valued shape parameter); see \link{fgev}
    }
    \item Distribution for random variables defined on bounded intervals:
    \itemize{
      \item \code{"uniform"} or \code{"unif"}: \code{min}, \code{max} (lower and upper boundaries), \code{lmass}, \code{umass} (point mass in lower or upper boundary); see \link{Uniform}
      \item \code{"beta"}: \code{shape1}, \code{shape2} (positive parameters); see \link{Beta}
    }
    \item Distributions for random variables with discrete / infinite support:
    \itemize{
      \item \code{"poisson"} or \code{"pois"}: \code{lambda} (positive mean); see \link{Poisson}
      \item \code{"negative-binomial"} or \code{"nbinom"}: \code{size} (positive dispersion parameter), \code{prob} (success probability), \code{mu} (mean, alternative to \code{prob}); see \link{NegBinomial}
    }
}
All numerical arguments should be of the same length. An exception are scalars of length 1, which will be recycled.
}
\examples{
crps(y = 1, family = "normal", mean = 0, sd = 2)
logs(y = 1, family = "normal", mean = 0, sd = 2)

crps(y = rnorm(20), family = "normal", mean = 1:20, sd = sqrt(1:20))
logs(y = rnorm(20), family = "normal", mean = 1:20, sd = sqrt(1:20))

## Arguments can have different lengths:
crps(y = rnorm(20), family = "normal", mean = 0, sd = 2)
crps(y = 1, family = "normal", mean = 1:20, sd = sqrt(1:20))

## Lists are accepted as input if elements with the required names are present:
df <- data.frame(y = rnorm(20), mean = 1:20, sd = 1)
crps(df, family = "normal")
logs(df, family = "normal")

## Mixture of normal distributions requires matrix input for parameters:
mval <- matrix(rnorm(20*50), nrow = 20)
sdval <- matrix(runif(20*50, min = 0, max = 2), nrow = 20)
weights <- matrix(rep(1/50, 20*50), nrow = 20)
crps(y = rnorm(20), family = "mixnorm", m = mval, s = sdval, w = weights)
}