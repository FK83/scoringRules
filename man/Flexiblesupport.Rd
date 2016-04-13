\name{Supplementary distributions: Flexible support}
\alias{fgev}
\alias{fgpd}
\alias{ftnorm}
\title{
Supplementary distributions (not in base R) with flexible support.
}
\description{We include the probability density functions of some distributions which are part of scoringRules, but are not part of base R. The parametrizations used here are identical to the ones used when calling \code{crps}, \code{logs} and \code{qs}.

Here we document distributions with flexible support: \code{fgev} - generalized extreme value distribution; \code{fgdp} - generalized Pareto distribution; \code{ftnorm} - truncated normal distribution.}
\usage{
fgev(x, location, scale, shape)
fgpd(x, location, scale, shape)
ftnorm(x, m, s, lb)
}
\arguments{
  \item{x}{vector of quantiles}
  \item{location}{vector of location parameters}
  \item{scale}{vector of scale parameters}
  \item{shape}{vector of shape parameters}
  \item{m}{vector of location parameters for \code{ftnorm}}
  \item{s}{vector of scale parameters for \code{ftnorm}}
  \item{lb}{lower truncation bound for \code{ftnorm}}
}
\details{
For details on generalized extreme value and generalized Pareto distributions, see Friederichs, F. and T.L. Thorarinsdottir (2012, "Forecast verification for extreme value distributions with an application to probabilistic peak wind prediction", Environmetrics 23, 579-594).

The truncated normal distribution is defined as the section of a normal distribution, rescaled to a probability distribution by an affine transformation (as opposed to a censored distribution which introduces point mass at the bounds).

The most common use in practice lies in the context of non-negative quantities, where a cut is made at the lower bound 0, leading to a pdf of
\deqn{f(x) = \phi(x)/(1 - \Phi(0)),} for \eqn{x \ge 0} and 0 otherwise. This example is for the standard normal distribution with pdf \eqn{\phi} and cdf \eqn{\Phi}, but also applies to location-scale transformed versions.

}
\keyword{distributions}
