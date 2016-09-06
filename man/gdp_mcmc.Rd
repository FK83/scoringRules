\name{GDP data}
\alias{gdp}
\alias{gdp_mcmc}
\docType{data}
\keyword{datasets}
\title{
Data and forecasts for US GDP growth
}
\description{
Historical data and forecast distributions for the growth rate of US gross domestic product (GDP). The forecasts are generated from a Bayesian Markov Switching model as described in Section 5 of Krueger et al (2016). 
}

\format{
\code{gdp} is a data frame which contains the real-time data set used in Section 5 of KLTG (2016), with the following columns:
\itemize{
\item \code{dt} - date in question (e.g., \code{"2013Q2"} for the second quarter of 2013)
\item \code{vint} - data vintage (i.e., the date at which the realization was recorded); same format as \code{dt}
\item \code{val} - value of the GDP growth rate
}

\code{gdp_mcmc} is a list, whereby each element is a data frame. \code{gdp_mcmc$forecasts} contains the simulated forecast distributions. There are 20 columns (corresponding to quarters 2008:Q1 to 2012:Q4) and 20,000 rows (corresponding to simulation draws). \code{gdp_mcmc$actuals} contains the actual observations. There are 20 columns (again corresponding to quarterly dates) and a single row. 
}

\source{
The GDP growth rate is computed from real-time data provided by the Federal Reserve Bank of Philadelphia, \url{http://www.phil.frb.org/research-and-data/real-time-center/} (series code ``ROUTPUT'', second-vintage data). The same data also enters the model which is used to generate the forecast distribution. \emph{Disclaimer: The provider of the raw data takes no responsibility for the accuracy of the data posted here. Furthermore, the raw data may be revised over time, and the website linked above should be consulted for the official, most recent version.}

The model from which the forecast draws are generated is described in Section 5 of KLTG (2016). Forecasts are one quarter ahead (that is, they are based on data until the previous quarter).

}
\details{
The realizations in \code{gdp_mcmc$actuals} are also contained in \code{gdp}, based on the second available vintage for each date. For example, \code{gdp_mcmc$actuals$X2008Q1} is the entry in \code{gdp} for which \code{dt == "2008Q1"} and \code{vint == "2008Q3"}.
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

Krueger, F., S. Lerch, T.L. Thorarinsdottir and T. Gneiting (2016), ``Probabilistic forecasting and comparative model assessment based on Markov Chain Monte Carlo output'', working paper, Heidelberg Institute for Theoretical Studies, available at \url{http://arxiv.org/abs/1608.06802}.

}
