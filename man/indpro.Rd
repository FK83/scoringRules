\name{indpro}
\alias{indpro}
\docType{data}
\title{
Example data: US industrial production 
}
\description{
Realizations and forecasts of the growth rate of quarterly US Industrial Production.
}

\format{
List; each element is a data frame. \code{dat} contains realizations and dates. \code{m_draws}, \code{v_draws} and \code{x_draws} contain the forecast mean, variance and (simulated) observation in each Gibbs sampler iteration;
please see the package vignette for details.
}

\source{Federal Reserve Bank of St. Louis (2014): Federal Reserve Economic Data, \url{http://research.stlouisfed.org/fred2/}.

\emph{Note: The underlying FRED series is called INDPRO. We compute quarterly growth rates via the following steps: 1) Download the monthly IP index from FRED. 2) Compute quarterly averages of the index. 3) Compute log growth
rates relative to previous quarter. 4) Multiply by 100.}
}

\examples{
\dontrun{

# Load data
data(indpro)

# Histogram of forecast draws for 2014Q2
hist(indpro$x_draws[,1], main = "Forecast draws for 2014:Q2", xlab = "Value")

# Add vertical line at realizing value
abline(v = subset(indpro$dat, dt == "2014Q2")$val, lwd = 3)
text(x = 9.5, y = 1500, "Realization")

}
}



\keyword{datasets}
