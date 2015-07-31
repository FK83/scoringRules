################################################################################

list_lsFunctions <- list()

################################################################################
### bounded interval

# uniform
ls.unif <- function(y, min, max) dunif(y, min, max, log=TRUE)
list_lsFunctions$'uniform' <- "ls.unif"

# beta
ls.beta <- function(y, shape1, shape2) dbeta(y, shape1, shape2, log=TRUE)
list_lsFunctions$'beta' <- "ls.beta"


################################################################################
### real line

# laplace
ls.lapl <- function(y, location, scale) log(flapl(y, location, scale))
list_lsFunctions$'laplace' <- "ls.lapl"

# logistic
ls.logis <- function(y, location, scale) dlogis(y, location, scale, log=TRUE)
list_lsFunctions$'logistic' <- "ls.logis"

# normal
ls.norm <- function(y, mean, sd) dnorm(y, mean, sd, log=TRUE)
list_lsFunctions$'normal' <- "ls.norm"

# mixture of normals
ls.mixn <- function(y, m, s, w) sapply(seq_along(y), function(i) lsmixnC(w[i, ], m[i, ], s[i, ], y[i]))
list_lsFunctions$'normal-mixture' <- "ls.mixn"

# two-piece-normal
ls.2pnorm <- function(y, m, s1, s2) log(f2pnorm(y, m, s1, s2))
list_lsFunctions$'two-piece-normal' <- "ls.2pnorm"

# t
ls.t <- function(y, df, location, scale) log(ft(y, df, location, scale))
list_lsFunctions$'t' <- "ls.t"


################################################################################
### non-negative

# exponential
ls.exp <- function(y, rate) dexp(y, rate, log=TRUE)
list_lsFunctions$'exponential' <- "ls.exp"

# gamma
ls.gamma <- function(y, shape, scale) dgamma(y, shape, scale = scale, log=TRUE)
list_lsFunctions$'gamma' <- "ls.gamma"

# log-laplace
ls.llapl <- function(y, locationlog, scalelog) log(fllapl(y, locationlog, scalelog))
list_lsFunctions$'log-laplace' <- "ls.llapl"

# log-logistic
ls.llogis <- function(y, locationlog, scalelog) log(fllogis(y, locationlog, scalelog))
list_lsFunctions$'log-logistic' <- "ls.llogis"

# log-normal
ls.lnorm <- function(y, meanlog, sdlog) dlnorm(y, meanlog, sdlog, log=TRUE)
list_lsFunctions$'log-normal' <- "ls.lnorm"

# truncated-normal
ls.tn <- function(y, m, s, lb) log(ftn(y, m, s, lb))
list_lsFunctions$'truncated-normal' <- "ls.tn"


################################################################################
### variable support

# generalized pareto distribution
ls.gpd <- function(y, location, scale, shape) log(fgpd(y, location, scale, shape))
list_lsFunctions$'gpd' <- "ls.gpd"

# generalized extreme value distribution
ls.gev <- function(y, location, scale, shape) log(fgev(y, location, scale, shape))
list_lsFunctions$'gev' <- "ls.gev"