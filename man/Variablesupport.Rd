\name{Supplementary distributions: Variable support}
\alias{fexp}
\alias{fgev}
\alias{fgpd}
\alias{fnorm}
\alias{ft}
\alias{flogis}
\keyword{distributions}
\author{Alexander Jordan}
\title{
Supplementary distributions (not in base R) with variable support.
}
\description{We include the probability density functions of some distributions which are part of scoringRules, but are not part of base R. The parametrizations used here are identical to the ones used when calling \code{crps} and \code{logs}.

Here we document distributions with variable support: \code{fexp} - location-scale exponential distribution with a point mass on the lower boundary; \code{fgdp} - generalized Pareto distribution with a point mass on the lower boundary; \code{fgev} - generalized extreme value distribution; \code{fnorm}, \code{flogis}, \code{ft} - (normal/logistic/Student's t)-distribution with flexible domain and point masses on the boundaries.}
\usage{
fexp(x, location, scale, mass = 0, log = FALSE)
fgpd(x, location, scale, shape, mass = 0, log = FALSE)

fgev(x, location, scale, shape)

fnorm(x, location, scale, lower = -Inf, upper = Inf, lmass = 0, umass = 0, log = FALSE)
ft(x, df, location, scale, lower = -Inf, upper = Inf, lmass = 0, umass = 0, log = FALSE)
flogis(x, location, scale, lower = -Inf, upper = Inf, lmass = 0, umass = 0, log = FALSE)
}
\arguments{
  \item{x}{vector of quantiles}
  \item{df}{vector of degrees of freedom parameters}
  \item{location}{vector of location parameters}
  \item{scale}{vector of scale parameters (positive)}
  \item{shape}{vector of shape parameters}
  \item{mass}{vector of point masses in \code{location}}
  \item{lower}{vector of lower bounds}
  \item{upper}{vector of upper bounds}
  \item{lmass}{vector of point masses in \code{lower}, or strings "trunc" / "cens"}
  \item{umass}{vector of point masses in \code{upper}, or strings "trunc" / "cens"}
  \item{log}{logical; if \code{TRUE}, the log of the density is returned}
}
\value{Density function of the relevant distribution, evaluated at \code{x}. NOTE: For distributions involving a point mass (e.g., when \code{lmass = "cens"} in \code{fnorm}), the density functions do not integrate to one.}
\details{
For details on generalized extreme value and generalized Pareto distributions, see Friederichs, F. and T.L. Thorarinsdottir (2012, "Forecast verification for extreme value distributions with an application to probabilistic peak wind prediction", Environmetrics 23, 579-594). Note that the support of both distributions depends on the input parameters; see \url{https://en.wikipedia.org/wiki/Generalized_extreme_value_distribution} and \url{https://en.wikipedia.org/wiki/Generalized_Pareto_distribution}.

Sometimes truncated or censored versions of the normal distribution are used to model variables with a restricted domain (e.g. precipitation). We allow the flexible specification of lower and upper boundaries and point masses in those boundaries. The truncated normal distribution assumes no point masses (i.e. redistributes the cut-off) and can be specified using the string "trunc" instead of a numerical probability. In contrast, the censored distribution introduces a point mass at the bound in the amount of the cut-off. Here, the string "cens" may be used for \code{lmass} or \code{umass}. The most common use in practice lies in the context of non-negative quantities. For example, a truncated standard normal distribution (left truncation at zero) has pdf \eqn{f(x) = \phi(x)/(1 - \Phi(0)),} for \eqn{x \ge 0} and 0 otherwise. A censored standard normal distribution (left censoring at zero) has point mass \eqn{\Phi(0)} at zero, and density \eqn{\phi(x)} for \eqn{x > 0}. 

The location-scale family based on Student's t-distribution (\code{ft}) has mean \eqn{location} for \eqn{df > 1} and variance \eqn{df/(df-2) * scale^2} for \eqn{df > 2}. Note that the \link{crps} exists only for \eqn{df > 1}. For details, see \url{https://en.wikipedia.org/wiki/Student's_t-distribution#Non-standardized_Student.27s_t-distribution}.
}
\seealso{The documentation for \link{crps} contains the full list of distributions supported by scoringRules (includes the ones documented here, as well as many others).}
