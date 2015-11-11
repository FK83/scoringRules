\name{sRTruncatednormal}
\alias{sRTruncatednormal}
\alias{ftnorm}
\title{
The Truncated Normal Distribution
}
\description{Density function for the truncated normal distribution.}
\usage{
ftnorm(x, m, s, lb)
}
\arguments{
  \item{x}{vector of quantiles}
  \item{m}{vector of location parameters}
  \item{s}{vector of scale parameters}
  \item{lb}{lower truncation bound}
}
\details{
The truncated normal distribution is defined as the section of a normal distribution, rescaled to a probability distribution by an affine transformation, as opposed to a censored distribution which achieves normalization by introducing point masses at the bounds.

The most common use in practice lies in the context of non-negative quantities, where a cut is made at the lower bound 0, leading to a pdf of
\deqn{f(x) = \phi(x)/(1 - \Phi(0)),} for \eqn{x \ge 0} and 0 otherwise. This example is for the standard normal distribution with pdf \eqn{\phi} and cdf \eqn{\Phi}, but also applies to location-scale transformed versions.

[Reference to Gneiting and Thorarinsdottir (2010)?]
}
