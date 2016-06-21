\name{gdp_mcmc}
\alias{gdp_mcmc}
\docType{data}
\title{
Example data: Forecasts and actual values of US GDP growth
}
\description{
Forecast distributions and actual values for the growth rate of US gross domestic product (GDP), 2008:Q1 to 2012:Q4. The forecasts are generated from a Bayesian Markov Switching model as described in Section 5 of Krueger et al (2016). 
}

\format{
List; each element is a data frame. \code{forecasts} contains the simulated forecast distributions. There are 20 columns (corresponding to quarters 2008:Q1 to 2012:Q4) and 20,000 rows (corresponding to simulation draws). \code{actuals} contains the actual observations. There are 20 columns (again corresponding to quarterly dates) and a single row. 
}

\source{
The GDP growth rate is computed from real-time data provided by the Federal Reserve Bank of Philadelphia, \url{http://www.phil.frb.org/research-and-data/real-time-center/} (series code ``ROUTPUT'', second-vintage data). The same data also enters the model which is used to generate the forecast distribution. \emph{Disclaimer: The provider of the raw data takes no responsibility for the accuracy of the data posted here. Furthermore, the raw data may be revised over time, and the website linked above should be consulted for the official, most recent version.}

The model from which the forecast draws are generated is described in Section 5 of Krueger et al (2016). Forecasts are one quarter ahead (that is, they are based on data until the previous quarter).

}

\examples{
\dontrun{

# Load data
data(gdp_mcmc)

# Histogram of forecast draws for 2012Q4
hist(gdp_mcmc$forecasts[, "X2012Q4"], main = "Forecast draws for 2012:Q4", xlab = "Value")

# Add vertical line at realizing value
abline(v = gdp_mcmc$actuals[, "X2012Q4"], lwd = 3)

}
}

\references{

Krueger, F., S. Lerch, T.L. Thorarinsdottir and T. Gneiting (2016), ``Probabilistic forecasting and comparative model assessment based on MCMC output'', working paper, Heidelberg Institute for Theoretical Studies, \emph{in preparation}.

}


\keyword{datasets}
