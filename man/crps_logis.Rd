\name{crps_logis}
\alias{crps_logis}
\alias{crps_clogis}
\alias{crps_tlogis}
\alias{gradcrps_logis}
\alias{gradcrps_clogis}
\alias{gradcrps_tlogis}
\alias{hesscrps_logis}
\alias{hesscrps_clogis}
\alias{hesscrps_tlogis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculating the CRPS for the logistic distribution
%%  ~~function to do ... ~~
}
\description{
These functions calculate the CRPS and its gradient and Hessian with respect
to the parameters of a location-scale transformation for the standard logistic
distribution. Furthermore, the censoring transformation and
the truncation transformation may be imposed on the
location-scale transformed logistic distribution.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
## CRPS functions
crps_logis(y, location = 0, scale = 1)
crps_clogis(y, location = 0, scale = 1, lower = -Inf, upper = Inf)
crps_tlogis(y, location = 0, scale = 1, lower = -Inf, upper = Inf)

## gradient (location, scale) functions
gradcrps_logis(y, location = 0, scale = 1)
gradcrps_clogis(y, location = 0, scale = 1, lower = -Inf, upper = Inf)
gradcrps_tlogis(y, location = 0, scale = 1, lower = -Inf, upper = Inf)

## Hessian (location, scale) functions
hesscrps_logis(y, location = 0, scale = 1)
hesscrps_clogis(y, location = 0, scale = 1, lower = -Inf, upper = Inf)
hesscrps_tlogis(y, location = 0, scale = 1, lower = -Inf, upper = Inf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{observation}
  \item{location}{location parameter}
  \item{scale}{scale parameter}
  \item{lower}{lower bound}
  \item{upper}{upper bound}
}
\details{something on logistic (standard/censored/truncated)
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
