\name{sRTwopiecenormal}
\alias{sRTwopiecenormal}
\alias{f2pnorm}
\title{
The Two-Piece Normal Distribution
}
\description{Density function for the two-piece normal distribution.}
\usage{
f2pnorm(x, location, scale1, scale2)
}
\arguments{
  \item{x}{vector of quantiles}
  \item{location}{vector of location parameters (mode of distribution)}
  \item{scale1, scale2}{vector of scale parameters (standard deviation to the left and right of the mode, must be positive)}
}
\details{
These functions are implemented for internal use.
For details, see Box A of Wallis (2004, "An Assessment of Bank of England and National Institute Inflation Forecast Uncertainties", National Institute Economic Review).
}
