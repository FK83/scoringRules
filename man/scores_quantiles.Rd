% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scores_quantiles.R
\name{scores_quantiles}
\alias{scores_quantiles}
\alias{qs_quantiles}
\alias{ints_quantiles}
\alias{qs_sample}
\alias{ints_sample}
\title{Quantile and interval scores}
\usage{
qs_quantiles(y, x, alpha)

ints_quantiles(y, x_lower, x_upper, target_coverage)

qs_sample(y, dat, alpha, w = NULL, type = 7, show_messages = TRUE)

ints_sample(y, dat, target_coverage, w = NULL, type = 7, show_messages = TRUE)
}
\arguments{
\item{y}{vector of observations}

\item{x}{vector of quantile predictions}

\item{alpha}{quantile level of interest}

\item{x_lower, x_upper}{vector of quantile predictions (lower and upper endpoints of prediction intervals)}

\item{target_coverage}{target (i.e., nominal) coverage level of prediction interval}

\item{dat}{vector or matrix (depending on \code{y}; see details)
of simulation draws from forecast distribution.}

\item{w}{vector of observation weights (currently ignored)}

\item{type}{integer between 1 and 9 that is passed on to stats function \link[stats]{quantile} (specifies algorithm for 
empirical quantile estimation; defaults to 7)}

\item{show_messages}{logical; display of messages (does not affect warnings and errors).}
}
\value{
A vector of score values. Smaller values indicate better forecasts. Note that 
the interval score refers to the central prediction interval at level \code{target_coverage}.
}
\description{
Compute quantile and interval scores, for  
given quantile predictions
}
\details{
For a vector \code{y} of length n, \code{dat} should be given as a matrix
with n rows. If \code{y} has length 1, then \code{dat} may be a vector.
}
\examples{
# Example 1: Illustrate that interval score is proportional to sum of two quantile scores
target_coverage <- .8
# corresponding quantile levels
alpha_1 <- .5*(1-target_coverage)
alpha_2 <- 1-.5*(1-target_coverage)
# compute interval score
ints_quantiles(y = 1, x_lower = qnorm(alpha_1), 
x_upper = qnorm(alpha_2), target_coverage = target_coverage)
# compute sum of quantile scores (scaled by 2/(1-target_coverage))
(2/(1-target_coverage))*(qs_quantiles(y = 1, x = qnorm(alpha_1), alpha = alpha_1) + 
qs_quantiles(y = 1, x = qnorm(alpha_2), alpha = alpha_2))

# Example 2: Compare exact to simulated quantile forecast from standard normal distribution
qs_quantiles(y = 1, x = qnorm(.1), alpha = .1)
qs_sample(y = 1, dat = rnorm(500), alpha = .1)
}
\references{
Quantile score

Koenker, R. and G. Bassett (1978): `Regression quantiles', Econometrica 46, 33-50. \doi{https://doi.org/10.2307/1913643}

Interval score

Gneiting, T. and A.E. Raftery (2007):
`Strictly proper scoring rules, prediction and estimation',
Journal of the American Statistical Association 102, 359-378. \doi{10.1198/016214506000001437}
}
\seealso{
The syntax of \code{\link{qs_sample}} and \code{\link{ints_sample}} is analogous to the functions documented on \code{\link{scores_sample_univ}}.
}
