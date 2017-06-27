\name{crps_norm}
\alias{crps_norm}
\alias{crps_cnorm}
\alias{crps_tnorm}
\alias{gradcrps_norm}
\alias{gradcrps_cnorm}
\alias{gradcrps_tnorm}
\alias{hesscrps_norm}
\alias{hesscrps_cnorm}
\alias{hesscrps_tnorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculating the CRPS for the normal distribution
%%  ~~function to do ... ~~
}
\description{
These functions calculate the CRPS and its gradient and Hessian with respect
to the parameters of a location-scale transformation for the standard normal
distribution. Furthermore, the censoring transformation and
the truncation transformation may be imposed on the
location-scale transformed normal distribution.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
## CRPS functions
crps_norm(y, mean = 0, sd = 1)
crps_cnorm(y, location = 0, scale = 1, lower = -Inf, upper = Inf)
crps_tnorm(y, location = 0, scale = 1, lower = -Inf, upper = Inf)

## gradient (location, scale) functions
gradcrps_norm(y, mean = 0, sd = 1)
gradcrps_cnorm(y, location = 0, scale = 1, lower = -Inf, upper = Inf)
gradcrps_tnorm(y, location = 0, scale = 1, lower = -Inf, upper = Inf)

## Hessian (location, scale) functions
hesscrps_norm(y, mean = 0, sd = 1)
hesscrps_cnorm(y, location = 0, scale = 1, lower = -Inf, upper = Inf)
hesscrps_tnorm(y, location = 0, scale = 1, lower = -Inf, upper = Inf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{observation}
  \item{mean}{mean parameter}
  \item{sd}{standard deviation parameter}
  \item{location}{location parameter}
  \item{scale}{scale parameter}
  \item{lower}{lower bound}
  \item{upper}{upper bound}
}
\details{something on normal (standard/censored/truncated)
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Vector for CRPS functions

Matrix with column names corresponding to the respective partial derivatives
for gradient and Hessian functions
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{
Alexander Jordan
%%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{others
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
