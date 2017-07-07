\name{crps_t}
\alias{crps_t}
\alias{crps_ct}
\alias{crps_tt}
\alias{gradcrps_t}
\alias{gradcrps_ct}
\alias{gradcrps_tt}
\alias{hesscrps_t}
\alias{hesscrps_ct}
\alias{hesscrps_tt}
\title{
Calculating the CRPS for Student's \eqn{t}-distribution
}
\description{
These functions calculate the CRPS and its gradient and Hessian with respect
to the parameters of a location-scale transformed Student's \eqn{t}-distribution.
Furthermore, the censoring transformation and
the truncation transformation may be introduced on top of the
location-scale transformed Student's \eqn{t}-distribution.
}
\usage{
## CRPS functions
crps_t(y, df, location = 0, scale = 1)
crps_ct(y, df, location = 0, scale = 1, lower = -Inf, upper = Inf)
crps_tt(y, df, location = 0, scale = 1, lower = -Inf, upper = Inf)

## gradient (location, scale) functions
gradcrps_t(y, df, location = 0, scale = 1)
gradcrps_ct(y, df, location = 0, scale = 1, lower = -Inf, upper = Inf)
gradcrps_tt(y, df, location = 0, scale = 1, lower = -Inf, upper = Inf)

## Hessian (location, scale) functions
hesscrps_t(y, df, location = 0, scale = 1)
hesscrps_ct(y, df, location = 0, scale = 1, lower = -Inf, upper = Inf)
hesscrps_tt(y, df, location = 0, scale = 1, lower = -Inf, upper = Inf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{observation}
  \item{df}{degrees of freedom parameter}
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
