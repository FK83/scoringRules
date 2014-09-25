\name{Simulation Study}
\alias{ergdist}
\alias{sim.fw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functions to replicate the simulation study in the paper by KLTG (2015)
}
\description{
Functions to simulate from the DGP used in the paper.
}
\usage{
ergdist(s, n)
sim.fw(f, s, n, ndraws)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{f, s, n}{Scalar DGP parameters for the Fox/West 2011 model}
  \item{ndraws}{Scalar (number of obserations to draw)
}
}
\value{
\item{par.ig}{Parameters of the inverse Gamma distribution for the Fox/West 2011 model}
\item{sim}{Vector of simulated draws from the Fox/West Model}
\item{pdf, cdf}{Probability density function and cumulative distribution function from the ergodic distribution}
\item{v, df}{Variance parameter and degrees of freedom parameter for ergodic t distribution. Note that the mean of this distribution is zero, and the variance is \eqn{\frac{v\times df}{df - 2}}}
}
\references{
E.B. Fox and M. West (2011), "Autoregressive Models for Variance Matrices: Stationary Inverse Wishart Processes", working paper, Duke University.

Krueger, F., S. Lerch, T. Thorarinsdottir and T. Gneiting (2015), ``Evaluating Simulation based Predictive Distributions'', working paper, Heidelberg Institute for Theoretical Studies, \emph{in preparation}.
}
\author{
Fabian Krueger
}
\note{
Much more detail needed here!
}

\examples{
\dontrun{
library(kltg)
set.seed(223)

# Parameters
f = 0.9
s = 4
n = 10
ndraws = 10000

# Simulate from the process used by KLTG
v = sim.fw(f, s, n, ndraws)
y1 = rnorm(ndraws, mean = 0, sd = sqrt(v$sim))

# Compute true ergodic dist
y2 = ergdist(s, n)

# Compare the two in a simple plot
d = density(y1)
plot(d, main = "black = simulated, red = true")
lines(x = d$x, y = y2$pdf(d$x), col = "red")

}


}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{simulation study}