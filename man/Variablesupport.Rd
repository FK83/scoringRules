\name{Supplementary distributions: Variable support}
\alias{fgev}
\alias{fgpd}
\alias{ftnorm}
\alias{fcnorm}
\alias{fctnorm}
\alias{ftcnorm}
\alias{fgnorm}
\title{
Supplementary distributions (not in base R) with variable support.
}
\description{We include the probability density functions of some distributions which are part of scoringRules, but are not part of base R. The parametrizations used here are identical to the ones used when calling \code{crps}, \code{logs} and \code{qs}.

Here we document distributions with variable support: \code{fgev} - generalized extreme value distribution; \code{fgdp} - generalized Pareto distribution; \code{ftnorm} - truncated normal distribution; \code{fcnorm} - censored normal distribution; \code{fctnorm} - normal distribution censored at left and truncated at right; \code{ftcnorm} - normal distribution truncated at left and censored at right; \code{fgnorm} - affine transformation of normal distribution.}
\usage{
fgev(x, location, scale, shape)
fgpd(x, location, scale, shape)

ftnorm(x, location, scale, lower, upper)
fcnorm(x, location, scale, lower, upper)
fctnorm(x, location, scale, lower, upper)
ftcnorm(x, location, scale, lower, upper)
fgnorm(x, location, scale, a, b, lower, upper)
}
\arguments{
  \item{x}{vector of quantiles}
  \item{location}{vector of location parameters}
  \item{scale}{vector of scale parameters (positive)}
  \item{shape}{vector of shape parameters}
  \item{a, b}{vector of parameters for the affine transformation of the normal CDF}
  \item{lower}{vector of lower bounds}
  \item{upper}{vector of upper bounds}
}
\value{Probability density function of the relevant distribution, evaluated at \code{x}.}
\details{
For details on generalized extreme value and generalized Pareto distributions, see Friederichs, F. and T.L. Thorarinsdottir (2012, "Forecast verification for extreme value distributions with an application to probabilistic peak wind prediction", Environmetrics 23, 579-594). Note that the support of both distributions depends on the input parameters; see \url{https://en.wikipedia.org/wiki/Generalized_extreme_value_distribution} and \url{https://en.wikipedia.org/wiki/Generalized_Pareto_distribution}.

The truncated normal distribution is defined as the section of a normal distribution, rescaled to integrate to one. In contrast, the censored distribution introduces point mass at the bounds. The most common use in practice lies in the context of non-negative quantities. For example, a truncated standard normal distribution (left truncation at zero) has pdf \eqn{f(x) = \phi(x)/(1 - \Phi(0)),} for \eqn{x \ge 0} and 0 otherwise. A censored standard normal distribution (left censoring at zero) has point mass \eqn{\Phi(0)} at zero, and density \eqn{\phi(x)} for \eqn{x > 0}. 

\code{ftcnorm} and \code{fctnorm} allow for flexible combinations of censoring and truncation.
}
\seealso{The documentation for \link{crps} contains the full list of distributions supported by scoringRules (includes the ones documented here, as well as many others).}
\keyword{distributions}
