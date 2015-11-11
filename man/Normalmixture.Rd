\name{sRNormalmixture}
\alias{sRNormalmixture}
\alias{fmixnorm}
\title{
The Mixture of Normal Distributions
}
\description{Density function for a mixture of normal distributions.}
\usage{
fmixnorm(x, m, s, w)
}
\arguments{
  \item{x}{vector of quantiles}
  \item{m}{matrix of means}
  \item{s}{matrix of standard deviations}
  \item{w}{matrix of weights}
  Each row corresponds to an entry in \code{x}. Each column corresponds to a mixture component. For convenience, scalars are allowed instead of matrices. Weights are normalized automatically.
}
\details{
The density function of a mixture of normal distributions is given by the weighted sum over the mixture components, \deqn{f(x) = \sum w_i/s_i \phi((x - m_i)/s_i),} where \eqn{\phi} is the pdf of the standard normal distribution.

''''
Matrices \code{m}, \code{s}, and \code{w}, where the \eqn{i}-th column of each matrix corresponds to the vector of mean values (\code{m[,i]}), standard deviations (\code{s[,i]}) and weights (\code{w[,i]}) of mixture component \eqn{i}. The weights in every row of \code{w} must sum to 1.
}
