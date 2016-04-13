\name{Supplementary distributions: Variable support}
\alias{fgev}
\alias{fgpd}
\alias{ftnorm}
\title{
Supplementary distributions (not in base R) with variable support.
}
\description{We include the probability density functions of some distributions which are part of scoringRules, but are not part of base R. The parametrizations used here are identical to the ones used when calling \code{crps}, \code{logs} and \code{qs}.

Here we document distributions with variable support: \code{fgev} - generalized extreme value distribution; \code{fgdp} - generalized Pareto distribution; \code{ftnorm} - truncated normal distribution.}
\usage{
fgev(x, location, scale, shape)
fgpd(x, location, scale, shape)
ftnorm(x, m, s, lb)
}
\arguments{
  \item{x}{vector of quantiles}
  \item{location}{vector of location parameters}
  \item{scale}{vector of scale parameters (positive)}
  \item{shape}{vector of shape parameters}
  \item{m}{vector of location parameters for \code{ftnorm}}
  \item{s}{vector of scale parameters for \code{ftnorm} (positive)}
  \item{lb}{lower truncation bound for \code{ftnorm}}
}
\value{Probability density function of the relevant distribution, evaluated at \code{x}.}
\details{
For details on generalized extreme value and generalized Pareto distributions, see Friederichs, F. and T.L. Thorarinsdottir (2012, "Forecast verification for extreme value distributions with an application to probabilistic peak wind prediction", Environmetrics 23, 579-594). Note that the support of both distributions depends on the input parameters; see \url{https://en.wikipedia.org/wiki/Generalized_extreme_value_distribution} and \url{https://en.wikipedia.org/wiki/Generalized_Pareto_distribution}.

The truncated normal distribution is defined as the section of a normal distribution, rescaled to a probability distribution by an affine transformation (as opposed to a censored distribution which introduces point mass at the bounds).

The most common use in practice lies in the context of non-negative quantities, where a cut is made at the lower bound 0, leading to a pdf of
\deqn{f(x) = \phi(x)/(1 - \Phi(0)),} for \eqn{x \ge 0} and 0 otherwise. This example is for the standard normal distribution with pdf \eqn{\phi} and cdf \eqn{\Phi}, but also applies to location-scale transformed versions.

}
\seealso{The documentation for \link{crps} contains the full list of distributions supported by scoringRules (includes the ones documented here, as well as many others).}
\keyword{distributions}
