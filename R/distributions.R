################################################################################
synonyms <- list(
  poisson = "pois",
  'negative-binomial' = "nbinom",
  uniform = "unif",
  laplace = "lapl",
  logistic = "logis",
  normal = "norm",
  'normal-mixture' = "mixnorm",
  'two-piece-normal' = "2pnorm",
  exponential = "exp",
  'log-laplace' = "llapl",
  'log-logistic' = "llogis",
  'log-normal' = "lnorm",
  'truncated-normal' = "tnorm",
  'censored-normal' = "cnorm",
  'truncated-censored-normal' = "tcnorm",
  'censored-truncated-normal' = "ctnorm"
)

################################################################################

# for only one choice of parameterization
checkNames1 <- function(input, reqinput) {
  ind <- match(reqinput, names(input), nomatch = 0)
  if (any(ind == 0)) {
    stop(paste("Missing parameter.",
               paste("Given input:", paste(names(input), collapse=", ")),
               paste("Required input:", paste(reqinput, collapse=", ")),
               sep="\n")
    )
  }
}

# for multiple parameterizations
checkNames2 <- function(input, reqinput) {
  ind <- lapply(reqinput, match, names(input), nomatch = 0)
  param.choice <- which(sapply(ind, function(x) all(x != 0)))
  if (length(param.choice) > 1) {
    stop(paste("Multiple parameterizations given. Please choose one.",
               paste("Given input:", paste(names(input), collapse=", ")),
               paste("Required input:", paste(lapply(reqinput, paste, collapse = ", "), collapse = "  OR  ")),
               sep="\n")
    ) 
  } else if (length(param.choice) == 0) {
    stop(paste("Missing parameter.",
               paste("Given input:", paste(names(input), collapse=", ")),
               paste("Required input:", paste(lapply(reqinput, paste, collapse = ", "), collapse = "  OR  ")),
               sep="\n")
    )
  } else {
    return(param.choice)
  }
}

# for scalar or vectorial parameters
checkVector <- function(input) {
  input_numeric <- sapply(input, is.numeric)
  if (any(!input_numeric)) {
    stop(paste("Non-numeric input:", paste(names(input)[!input_numeric], collapse=", ")))
  }
  input_isvector <- sapply(input, is.vector)
  if (any(!input_isvector)) {
    stop(paste("Non-scalar or non-vectorial input:", paste(names(input)[!input_isvector], collapse=", ")))
  }
  input_lengths <- sapply(input, length)
  input_maxlength <- max(input_lengths)
  if (any(is.na(match(input_lengths, c(1, input_maxlength))))) {
    stop("Incompatible input vector lengths. Each vector should be of length 1 or n.")
  }
}

# for matrix parameters
checkMatrix <- function(input) {
  input_numeric <- sapply(input, is.numeric)
  if (any(!input_numeric)) {
    stop(paste("Non-numeric input:", paste(names(input)[!input_numeric], collapse=", ")))
  }
  if (!is.vector(input$y)) stop("Non-scalar or non-vectorial observation 'y' input.")
  input_ismatrix <- sapply(input[-1], is.matrix)
  if (any(!input_ismatrix)) {
    stop(paste("Non-matrix input:", paste(names(input[-1])[!input_ismatrix], collapse=", ")))
  }
  input_dims <- sapply(input[-1], dim)
  ylength <- length(input$y)
  input_maxdim <- c(ylength, max(input_dims[2,]))
  if (any(is.na(sapply(1:2, function(i) match(input_dims[i,], input_maxdim[i]))))) {
    stop("Incompatible input object sizes.")
  }
}

################################################################################
### discrete / infinite support

### poisson
check.pois <- function(input) {
  reqinput <- c("y", "lambda")
  checkNames1(input, reqinput)
  input <- input[reqinput]
  checkVector(input)
  
  if (any(!input$lambda>0)) stop("Parameter 'lambda' contains non-positive values.")
  if (any(is.infinite(input$lambda))) stop("Parameter 'lambda' contains infinite values.")
  
  return(input)
}

### negative binomial
check.nbinom <- function(input) {
  reqinput <- list(c("y", "size", "prob"),
                   c("y", "size", "mu"))
  choice <- checkNames2(input, reqinput)
  input <- input[reqinput[[choice]]]
  checkVector(input)
  
  if (choice == 1) {
    if (any(input$prob > 1 | input$prob <= 0)) stop("Parameter 'prob' not in (0, 1]")
  } else if (choice == 2) {
    if (any(input$mu < 0)) stop("Parameter 'mu' contains negative values.")
    if (any(is.infinite(input$mu))) stop("Parameter 'mu' contains infinite values.")
    input$prob <- input$size/(input$size + input$mu)
  }
  if (any(!input$size>0)) stop("Parameter 'size' contains non-positive values.")
  if (any(!is.finite(input$size))) stop("Parameter 'size' contains infinite values.")
  
  return(input[c("y", "size", "prob")])
}

################################################################################
### bounded interval

### uniform
check.unif <- function(input) {
  reqinput <- c("y", "min", "max")
  checkNames1(input, reqinput)
  input <- input[reqinput]
  checkVector(input)
  
  if (any(is.infinite(c(input$min, input$max)))) stop("Invalid distribution due to infinite bounds.")
  if (any(input$min > input$max)) stop("Parameter 'min' contains greater values than parameter 'max'.")
  
  return(input)
}

### beta
check.beta <- function(input) {
  reqinput <- c("y", "shape1", "shape2")
  checkNames1(input, reqinput)
  input <- input[reqinput]
  checkVector(input)
  
  if (any(!input$shape1>0)) stop("Parameter 'shape1' contains non-positive values.")
  if (any(is.infinite(input$shape1))) stop("Parameter 'shape1' contains infinite values.")
  if (any(!input$shape2>0)) stop("Parameter 'shape2' contains non-positive values.")
  if (any(is.infinite(input$shape2))) stop("Parameter 'shape2' contains infinite values.")
  
  return(input)
}

################################################################################
### real line

### laplace
check.lapl <- function(input) {
  reqinput <- c("y", "location", "scale")
  checkNames1(input, reqinput)
  input <- input[reqinput]
  checkVector(input)
  
  if (any(is.infinite(input$location))) stop("Parameter 'location' contains infinite values.")
  if (any(!input$scale>0)) stop("Parameter 'scale' contains non-positive values.")
  if (any(is.infinite(input$scale))) stop("Parameter 'scale' contains infinite values.")
  
  return(input)
}

flapl <- function(x, location, scale) {
  dexp(abs(x - location), 1/scale) / 2
}
f2plapl <- function(x, location, scale1, scale2) {
  n <- max(length(x), length(location), length(scale1), length(scale2))
  z <- rep(x - location, len = n)
  s <- ifelse(z < 0, scale1, scale2)
  s / (scale1 + scale2) * dexp(abs(z), 1/s)
}

### logistic
check.logis <- function(input) {
  reqinput <- c("y", "location", "scale")
  checkNames1(input, reqinput)
  input <- input[reqinput]
  checkVector(input)
  
  if (any(is.infinite(input$location))) stop("Parameter 'location' contains infinite values.")
  if (any(!input$scale>0)) stop("Parameter 'scale' contains non-positive values.")
  if (any(is.infinite(input$scale))) stop("Parameter 'scale' contains infinite values.")
  
  return(input)
}

### normal
check.norm <- function(input) {
  reqinput <- list(c("y", "mean", "sd"),
                   c("y", "location", "scale"))
  choice <- checkNames2(input, reqinput)
  input <- input[reqinput[[choice]]]
  checkVector(input)
  
  if (choice == 2) {
    names(input)[names(input) == "location"] <- "mean"
    names(input)[names(input) == "scale"] <- "sd"
    if (any(is.infinite(input$mean))) stop("Parameter 'location' contains infinite values.")
    if (any(!input$sd>0)) stop("Parameter 'scale' contains non-positive values.")
    if (any(is.infinite(input$sd))) stop("Parameter 'scale' contains infinite values.")
  } else if (choice == 1) {
    if (any(is.infinite(input$mean))) stop("Parameter 'mean' contains infinite values.")
    if (any(!input$sd>0)) stop("Parameter 'sd' contains non-positive values.")
    if (any(is.infinite(input$sd))) stop("Parameter 'sd' contains infinite values.")
  }
  
  return(input)
}

### normal-mixture
check.mixnorm <- function(input) {
  reqinput <- c("y", "m", "s", "w")
  checkNames1(input, reqinput)
  input <- input[reqinput]
  checkMatrix(input)
  
  if (any(is.infinite(input$m))) stop("Parameter 'm' contains infinite values.")
  if (any(!input$s>0)) stop("Parameter 's' contains non-positive values.")
  if (any(is.infinite(input$s))) stop("Parameter 's' contains infinite values.")
  if (any(input$w < 0 | input$w > 1)) stop("Parameter 'w' contains values not in [0, 1].")
  if (all.equal(apply(input$w, 1, sum), rep(1, dim(input$w)[1])) != TRUE) stop("Parameter 'w' contains weighting schemes which do not sum up to 1.")
  
  return(input)
}
fmixnorm <- function(x, m, s, w) {
  if (!is.vector(x)) stop("object 'x' not a vector")
  p <- length(x)
  
  if (is.matrix(m)) {
    m.q <- dim(m)[2]
    if (dim(m)[1] != p) stop("number of rows of object 'm' does not match length of object 'x'")
  } else if (length(m) != 1) {
    stop("object 'm' is neither a matrix nor a single number")
  }
  
  if (is.matrix(s)) {
    s.q <- dim(s)[2]
    if (dim(s)[1] != p) stop("number of rows of object 's' does not match length of object 'x'")
    if (exists("m.q")) {
      if (s.q != m.q) stop("number of mixture components in 's' and 'm' do not match")
    }
  } else if (length(s) != 1) {
    stop("object 's' is neither a matrix nor a single number")
  }
  
  if (is.matrix(w)) {
    w.q <- dim(w)[2]
    if (dim(w)[1] != p) stop("number of rows of object 'w' does not match length of object 'x'")
    if (exists("m.q")) {
      if (w.q != m.q) stop("number of mixture components in 'w' and 'm' do not match")
    }
    if (exists("s.q")) {
      if (w.q != s.q) stop("number of mixture components in 'w' and 's' do not match")
    }
  } else if (length(w) == 1) {
    w <- 1/max(dim(m)[2], dim(s)[2], 1)
  } else {
    stop("object 'w' is neither a matrix nor a single number")
  }
  
  rowSums(dnorm(as.matrix(x), m, s) * w)
}

### two-piece-normal
check.2pnorm <- function(input) {
  reqinput <- c("y", "location", "scale1", "scale2")
  checkNames1(input, reqinput)
  input <- input[reqinput]
  checkVector(input)
  
  if (any(is.infinite(input$location))) stop("Parameter 'location' contains infinite values.")
  if (any(!input$scale1 > 0)) stop("Parameter 'scale1' contains non-positive values.")
  if (any(is.infinite(input$scale1))) stop("Parameter 'scale1' contains infinite values.")
  if (any(!input$scale2>0)) stop("Parameter 'scale2' contains non-positive values.")
  if (any(is.infinite(input$scale2))) stop("Parameter 'scale2' contains infinite values.")
  
  return(input)
}
f2pnorm <- function(x, location, scale1, scale2) {
  n <- max(length(x), length(location), length(scale1), length(scale2))
  z <- rep(x - location, len = n)
  s <- ifelse(z < 0, scale1, scale2)
  2 * s / (scale1 + scale2) * dnorm(z, 0, s)
}

### t
check.t <- function(input) {
  reqinput <- c("y", "df", "location", "scale")
  checkNames1(input, reqinput)
  input <- input[reqinput]
  checkVector(input)
  
  if (any(!input$df>0)) stop("Parameter 'df' contains non-positive values.")
  if (any(is.infinite(input$location))) stop("Parameter 'location' contains infinite values.")
  if (any(!input$scale>0)) stop("Parameter 'scale' contains non-positive values.")
  if (any(is.infinite(input$scale))) stop("Parameter 'scale' contains infinite values.")
  
  return(input)
}
ft <- function(x, df, location, scale) {
  z <- (x - location) / scale
  1/scale * dt(z, df)
}


################################################################################
### non-negative

### exponential
check.exp <- function(input) {
  reqinput <- c("y", "rate")
  checkNames1(input, reqinput)
  input <- input[reqinput]
  checkVector(input)
  
  if (any(!input$rate>0)) stop("Parameter 'rate' contains non-positive values.")
  if (any(is.infinite(input$rate))) stop("Parameter 'rate' contains infinite values.")
  
  return(input)
}

### gamma
check.gamma <- function(input) {
  reqinput <- list(c("y", "shape", "rate"),
                   c("y", "shape", "scale")
  )
  choice <- checkNames2(input, reqinput)
  input <- input[reqinput[[choice]]]
  checkVector(input)
  
  if (any(!input$shape>0)) stop("Parameter 'shape' contains non-positive values.")
  if (any(!is.finite(input$shape))) stop("Parameter 'shape' contains infinite values.")
  if (choice == 1) {
    if (any(!input$rate>0)) stop("Parameter 'rate' contains non-positive values.")
    if (any(!is.finite(input$rate))) stop("Parameter 'rate' contains infinite values.")
    input$scale <- 1/input$rate
  } else if (choice == 2) {
    if (any(!input$scale>0)) stop("Parameter 'scale' contains non-positive values.")
    if (any(!is.finite(input$scale))) stop("Parameter 'scale' contains infinite values.")
  }
  
  return(input[c("y", "shape", "scale")])
}

### log-laplace
check.llapl <- function(input) {
  reqinput <- c("y", "locationlog", "scalelog")
  checkNames1(input, reqinput)
  input <- input[reqinput]
  checkVector(input)
  
  if (any(is.infinite(input$locationlog))) stop("Parameter 'locationlog' contains infinite values.")
  if (any(!input$scalelog>0)) stop("Parameter 'scalelog' contains non-positive values.")
  if (any(is.infinite(input$scalelog))) stop("Parameter 'scalelog' contains infinite values.")
  
  return(input)
}
fllapl <- function(x, locationlog, scalelog) {
  x1 <- log(pmax(x, 0))
  ind <- is.infinite(x1)
  d <- 1/x * flapl(x1, locationlog, scalelog)
  d[ind] <- 0
  return(d)
}

### log-logistic
check.llogis <- function(input) {
  reqinput <- c("y", "locationlog", "scalelog")
  checkNames1(input, reqinput)
  input <- input[reqinput]
  checkVector(input)
  
  if (any(is.infinite(input$locationlog))) stop("Parameter 'locationlog' contains infinite values.")
  if (any(!input$scalelog>0)) stop("Parameter 'scalelog' contains non-positive values.")
  if (any(is.infinite(input$scalelog))) stop("Parameter 'scalelog' contains infinite values.")
  
  return(input)
}
fllogis <- function(x, locationlog, scalelog) {
  x1 <- log(pmax(x, 0))
  ind <- is.infinite(x1)
  d <- 1/x * dlogis(x1, locationlog, scalelog)
  d[ind] <- 0
  return(d)
}

### log-normal
check.lnorm <- function(input) {
  reqinput <- list(c("y", "meanlog", "sdlog"),
                   c("y", "locationlog", "scalelog"))
  choice <- checkNames2(input, reqinput)
  input <- input[reqinput[[choice]]]
  checkVector(input)
  
  if (choice == 2) {
    names(input)[names(input) == "locationlog"] <- "meanlog"
    names(input)[names(input) == "scalelog"] <- "sdlog"
    if (any(is.infinite(input$meanlog))) stop("Parameter 'locationlog' contains infinite values.")
    if (any(!input$sdlog>0)) stop("Parameter 'scalelog' contains non-positive values.")
    if (any(is.infinite(input$sdlog))) stop("Parameter 'scalelog' contains infinite values.")
  } else {
    if (any(is.infinite(input$meanlog))) stop("Parameter 'meanlog' contains infinite values.")
    if (any(!input$sdlog>0)) stop("Parameter 'sdlog' contains non-positive values.")
    if (any(is.infinite(input$sdlog))) stop("Parameter 'sdlog' contains infinite values.")
  }
  
  return(input)
}

### truncated-normal
check.tnorm <- function(input) {
  reqinput <- c("y", "location", "scale")
  optinput <- c("lower", "upper")
  checkNames1(input, reqinput)
  input <- input[c(reqinput, optinput[optinput %in% names(input)])]
  checkVector(input)
  
  if (any(is.infinite(input$location))) stop("Parameter 'location' contains infinite values.")
  if (any(!input$scale > 0)) stop("Parameter 'scale' contains non-positive values.")
  if (any(is.infinite(input$scale))) stop("Parameter 'scale' contains infinite values.")
  
  return(input)
}
ftnorm <- function(x, location, scale, lower, upper) {
  a <- 1 - pnorm(upper, location, scale, lower.tail = FALSE) - pnorm(lower, location, scale)
  d <- dnorm(x, location, scale) / a
  d[x < lower] <- 0
  d[x > upper] <- 0
  return(d)
}

### censored-normal
check.cnorm <- function(input) {
  check.tnorm(input)
}
fcnorm <- function(x, location, scale, lower, upper) {
  d <- dnorm(x, location, scale)
  d[x < lower] <- 0
  d[x > upper] <- 0
  return(d)
}

### truncated-censored normal
check.ctnorm <- function(input) {
  check.tnorm(input)
}
fctnorm <- function(x, location, scale, lower, upper) {
  a <- 1 - pnorm(upper, location, scale, lower.tail = FALSE)
  d <- dnorm(x, location, scale) / a
  d[x < lower] <- 0
  d[x > upper] <- 0
  return(d)
}

### truncated-censored normal
check.tcnorm <- function(input) {
  check.tnorm(input)
}
ftcnorm <- function(x, location, scale, lower, upper) {
  a <- 1 - pnorm(lower, location, scale)
  d <- dnorm(x, location, scale) / a
  d[x < lower] <- 0
  d[x > upper] <- 0
  return(d)
}

check.gnorm <- function(input) {
  reqinput <- c("y", "location", "scale", "a", "b", "lower", "upper")
  checkNames1(input, reqinput)
  input <- input[reqinput]
  checkVector(input)
  
  if (any(is.infinite(input$location))) stop("Parameter 'location' contains infinite values.")
  if (any(!input$scale > 0)) stop("Parameter 'scale' contains non-positive values.")
  if (any(is.infinite(input$scale))) stop("Parameter 'scale' contains infinite values.")
  
  return(input)
}
fgnorm <- function(x, location, scale, a, b, lower, upper) {
  d <- dnorm(x, location, scale) * a
  d[x < lower] <- 0
  d[x > upper] <- 0
  return(d)
}

################################################################################
### variable support

### gpd
check.gpd <- function(input) {
  reqinput <- c("y", "location", "scale", "shape")
  checkNames1(input, reqinput)
  input <- input[reqinput]
  checkVector(input)
  
  if (any(is.infinite(input$location))) stop("Parameter 'location' contains infinite values.")
  if (any(!input$scale>0)) stop("Parameter 'scale' contains non-positive values.")
  if (any(is.infinite(input$scale))) stop("Parameter 'scale' contains infinite values.")
  if (any(is.infinite(input$shape))) stop("Parameter 'shape' contains infinite values.")
  
  return(input)
}
fgpd <- function(x, location, scale, shape) {
  z <- (x - location) / scale
  
  ind <- abs(shape) < 1e-12
  if (any(ind)) {
    if (length(z) < length(shape))
      z <- rep(z, len = length(shape))
    if (length(scale) < length(z))
      scale <- rep(scale, len = length(z))
    out <- numeric(length(z))
    out[ind] <- 1/scale[ind] * dexp(z[ind], 1)
    out[!ind] <- 1/scale[!ind] * fgpd(z[!ind], 0, 1, shape[!ind])
  } else {
    out <- 1/scale * (1 + shape * z)^(- 1 - 1/shape)
    out[z < 0] <- 0
    out[z > -1/shape & shape < 0] <- 0
  }
  
  return(out)
}

### gev
check.gev <- function(input) {
  reqinput <- c("y", "location", "scale", "shape")
  checkNames1(input, reqinput)
  input <- input[reqinput]
  checkVector(input)
  
  if (any(is.infinite(input$location))) stop("Parameter 'location' contains infinite values.")
  if (any(!input$scale>0)) stop("Parameter 'scale' contains non-positive values.")
  if (any(is.infinite(input$scale))) stop("Parameter 'scale' contains infinite values.")
  if (any(is.infinite(input$shape))) stop("Parameter 'shape' contains infinite values.")
      
  return(input)
}

fgev <- function(x, location, scale, shape) {
  z <- (x - location) / scale
  
  ind <- abs(shape) < 1e-12
  if (any(ind)) {
    if (length(z) < length(shape))
      z <- rep(z, len = length(shape))
    if (length(scale) < length(z))
      scale <- rep(scale, len = length(z))
    out <- numeric(length(z))
    out[ind] <- 1 / scale[ind] * exp(-z[ind] - exp(-z[ind]))
    out[!ind] <- 1 / scale[!ind] * fgev(z[!ind], 0, 1, shape[!ind])
  } else {
    zz <- 1 + shape * z
    out <-
      1 / scale * exp(log1p(shape * z) * (-1 - 1 / shape) - zz ^ (-1 / shape))
    out[z * sign(shape) < -1 / abs(shape)] <- 0
  }
  
  return(out)
}
