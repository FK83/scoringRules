\name{crps_norm}
\alias{crps_norm}
\alias{crps_cnorm}
\alias{crps_tnorm}
\alias{gradcrps_norm}
\alias{gradcrps_cnorm}
\alias{gradcrps_tnorm}
\alias{hesscrps_norm}
\alias{hesscrps_cnorm}
\alias{hesscrps_tnorm}
\title{
Calculating the CRPS for the normal distribution
}
\description{
These functions calculate the CRPS and its gradient and Hessian with respect
to the parameters of a location-scale transformed normal
distribution. Furthermore, the censoring transformation and
the truncation transformation may be introduced on top of the
location-scale transformed normal distribution.
}
\usage{
## CRPS functions
crps_norm(y, location = 0, scale = 1)
crps_cnorm(y, location = 0, scale = 1, lower = -Inf, upper = Inf)
crps_tnorm(y, location = 0, scale = 1, lower = -Inf, upper = Inf)

## gradient (location, scale) functions
gradcrps_norm(y, location = 0, scale = 1)
gradcrps_cnorm(y, location = 0, scale = 1, lower = -Inf, upper = Inf)
gradcrps_tnorm(y, location = 0, scale = 1, lower = -Inf, upper = Inf)

## Hessian (location, scale) functions
hesscrps_norm(y, location = 0, scale = 1)
hesscrps_cnorm(y, location = 0, scale = 1, lower = -Inf, upper = Inf)
hesscrps_tnorm(y, location = 0, scale = 1, lower = -Inf, upper = Inf)
}
\arguments{
  \item{y}{observation}
  \item{location}{location parameter}
  \item{scale}{scale parameter}
  \item{lower}{lower bound}
  \item{upper}{upper bound}
}
\value{
For the CRPS functions: a vector of score values.

For the gradient and Hessian functions: a matrix with column names corresponding to the respective partial derivatives.
}
\author{
Alexander Jordan
}
