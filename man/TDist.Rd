\name{sRTDist}
\alias{sRTDist}
\alias{ft}
\title{
The Student t Distribution (location-scale family)
}
\description{Density function for the location-scale family based on Student's t-distribution.}
\usage{
ft(x, df, location, scale)
}
\arguments{
  \item{x}{vector of quantiles}
  \item{df}{vector of degrees of freedom}
  \item{location}{vector of location parameters}
  \item{scale}{vector of scale parameters}
}
\details{
The location-scale family based on Student's t-distribution has mean \eqn{location} for \eqn{df > 1} and variance \eqn{df/(df-2) * scale^2} for \eqn{df > 2}. Consequently, the \link{crps} exists only for \eqn{df > 1}.
}
