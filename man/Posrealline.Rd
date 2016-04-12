\name{Additional distributions on the positive real line}
\alias{fllapl}
\alias{fllogis}
\title{
Additional distributions (not in base R) supported on the positive real line.
}
\description{We include the probability density functions of some distributions which are part of scoringRules, but are not part of base R. The parametrizations used here are identical to the ones used when calling \code{crps}, \code{logs} and \code{qs}.

Here we document distributions on the positive real line: \code{fllapl} - log Laplace distribution; \code{fllogis} - log logistic distribution.}
\usage{
fllapl(x, locationlog, scalelog)
fllogis(x, locationlog, scalelog)
}
\arguments{
  \item{x}{vector of quantiles}
  \item{locationlog}{vector of location parameters on the log scale}
  \item{scalelog}{vector of scale parameters on the log scale}
}
\value{Probability density function of the relevant distribution, evaluated at \code{x}.}
\details{
To be added.
}
\seealso{The documentation for \link{crps} contains the full list of distributions supported by scoringRules (includes the ones documented here, as well as many others).}