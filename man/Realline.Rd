\name{Supplementary distributions: Real line}
\alias{flapl}
\alias{f2pexp}
\alias{fmixnorm}
\alias{f2pnorm}
\keyword{distributions}
\author{Alexander Jordan}
\title{Supplementary distributions (not in base R) supported on the real line.}
\description{We include the probability density functions of some distributions which are part of scoringRules, but are not part of base R. The parametrizations used here are identical to the ones used when calling \code{crps} and \code{logs}.

Here we document distributions with support on the real line: \code{flapl} - Laplace distribution; \code{f2pexp} - two-piece exponential distribution; \code{fmixnorm} - mixture of normal distributions; \code{f2pnorm} - two-piece normal distribution.
}
\usage{
flapl(x, location, scale)

f2pexp(x, location, scale1, scale2)
f2pnorm(x, location, scale1, scale2)

fmixnorm(x, m, s, w)
}
\arguments{
  \item{x}{vector of quantiles}
  \item{location}{vector of location parameters}
  \item{scale, scale1, scale2}{vector of scale parameters}
  \item{m}{matrix of means (rows correspond to observations, columns correspond to mixture components)}
  \item{s}{matrix of standard deviations (same structure as \code{m})}
  \item{w}{matrix of weights (same structure as \code{m})}
}
\value{Probability density function of the relevant distribution, evaluated at \code{x}.}
\details{
The Laplace distribution (\code{flapl}) is described on \url{https://en.wikipedia.org/wiki/Laplace_distribution}. It is a special case of the two-piece exponential distribution (\code{f2pexp}), which allows for different scale parameters to the left and right of \code{location}.

The density function of a mixture of normal distributions (\code{fmixnorm}) is given by the weighted sum over the mixture components, \deqn{f(x) = \sum w_i/s_i \phi((x - m_i)/s_i),} where \eqn{\phi} is the pdf of the standard normal distribution.

For details on the two-piece normal distribution (\code{f2pnorm}), see Box A of Wallis (2004, "An Assessment of Bank of England and National Institute Inflation Forecast Uncertainties", National Institute Economic Review).
}
\seealso{The documentation for \link{crps.numeric} contains the full list of distributions supported by scoringRules (includes the ones documented here, as well as many others).

\link{fnorm}, \link{flogis}, \link{ft}
}
\examples{
# Plot PDF of Laplace distribution
ff <- function(x) flapl(x, location = 0, scale = 2)
curve(ff, from = -8, to = 8, bty = "n", xlab = "Value", 
      ylab = "PDF", 
      main = "Laplace distribution with location 0 and scale 2")
}
