\name{crps_t}
\alias{crps_t}
\alias{crps_ct}
\alias{crps_tt}
\alias{gradcrps_t}
\alias{gradcrps_ct}
\alias{gradcrps_tt}
\alias{hesscrps_t}
\alias{hesscrps_ct}
\alias{hesscrps_tt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculating the CRPS for the t distribution
%%  ~~function to do ... ~~
}
\description{
These functions calculate the CRPS and its gradient and Hessian with respect
to the parameters of a location-scale transformation for the standard t
distribution. Furthermore, the censoring transformation and
the truncation transformation may be imposed on the
location-scale transformed t distribution.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
## CRPS functions
crps_t(y, df, location = 0, scale = 1)
crps_ct(y, df, location = 0, scale = 1, lower = -Inf, upper = Inf)
crps_tt(y, df, location = 0, scale = 1, lower = -Inf, upper = Inf)

## gradient (location, scale) functions
gradcrps_t(y, df, location = 0, scale = 1)
gradcrps_ct(y, df, location = 0, scale = 1, lower = -Inf, upper = Inf)
gradcrps_tt(y, df, location = 0, scale = 1, lower = -Inf, upper = Inf)

## Hessian (location, scale) functions
hesscrps_t(y, df, location = 0, scale = 1)
hesscrps_ct(y, df, location = 0, scale = 1, lower = -Inf, upper = Inf)
hesscrps_tt(y, df, location = 0, scale = 1, lower = -Inf, upper = Inf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{observation}
  \item{df}{degrees of freedom parameter}
  \item{location}{location parameter}
  \item{scale}{scale parameter}
  \item{lower}{lower bound}
  \item{upper}{upper bound}
}
\details{something on t (standard/censored/truncated)
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
