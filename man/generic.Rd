\name{generic}
\alias{crps}
\alias{logs}
\title{
Generic Scoring Rule Calculation
}
\description{
Generic functions for calculating the Continuous Ranked Probability Score and the Logarithmic Score of R objects.

\code{scoringRules} provides default methods (\link{crps.numeric}, \link{logs.numeric}) to calculate scores of forecasts that are members of families of parametric distributions.
}
\usage{
crps(y, ...)
logs(y, ...)
}
\arguments{
  \item{y}{an object for which the score is to be calculated}
  \item{...}{further arguments passed to or from other methods}
}
\value{
Returns a vector of scores. One for each forecast-observation pair.
}
\author{
Alexander Jordan, Fabian Krueger, Sebastian Lerch
}
\seealso{
\link{crps.numeric}, \link{logs.numeric}
}