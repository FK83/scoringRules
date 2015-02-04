################################################################################
list_inputChecks <- list()

# for no choice of parameterization
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
               paste("Required input:", paste(names(reqinput), collapse=", ")),
               sep="\n")
    ) 
  } else if (length(param.choice) == 0) {
    stop(paste("Missing parameter.",
               paste("Given input:", paste(names(input), collapse=", ")),
               paste("Required input:", paste(names(reqinput), collapse=", ")),
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
  if (anyNA(match(input_lengths, c(1, input_maxlength)))) {
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
  if (anyNA(sapply(1:2, function(i) match(input_dims[i,], input_maxdim[i])))) {
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
  
  for (i in seq_along(input)) {
    assign(names(input)[i], input[[i]])
  }
  if (any(!lambda>0)) stop("Parameter 'lambda' contains non-positive values.")
  if (any(is.infinite(lambda))) stop("Parameter 'lambda' contains infinite values.")
  
  return(list(y = y, lambda = lambda))
}
list_inputChecks$'poisson' <- "check.pois"

### negative binomial
check.nbinom <- function(input) {
  reqinput <- list(c("y", "size", "prob"),
                   c("y", "size", "mu"))
  choice <- checkNames2(input, reqinput)
  input <- input[reqinput[[choice]]]
  checkVector(input)
  
  for (i in seq_along(input)) {
    assign(names(input)[i], input[[i]])
  }
  if (choice == 1) {
    if (any(prob > 1 | prob <= 0)) stop("Parameter 'prob' not in (0, 1]")
  } else if (choice == 2) {
    if (any(mu < 0)) stop("Parameter 'mu' contains negative values.")
    if (any(is.infinite(mu))) stop("Parameter 'mu' contains infinite values.")
    prob <- size/(size + mu)
  }
  if (any(!size>0)) stop("Parameter 'size' contains non-positive values.")
  if (any(!is.finite(size))) stop("Parameter 'size' contains infinite values.")
  
  return(list(y = y, size = size, prob = prob))
}
list_inputChecks$'negative-binomial' <- "check.nbinom"

################################################################################
### bounded interval

### uniform
check.unif <- function(input) {
  reqinput <- c("y", "min", "max")
  checkNames1(input, reqinput)
  input <- input[reqinput]
  checkVector(input)
  
  for (i in seq_along(input)) {
    assign(names(input)[i], input[[i]])
  }
  if (any(is.infinite(c(min, max)))) stop("Invalid distribution due to infinite bounds.")
  if (any(min > max)) stop("Parameter 'min' contains greater values than parameter 'max'.")
  
  return(list(y = y, min = min, max = max))
}
list_inputChecks$'uniform' <- "check.unif"

### beta
check.beta <- function(input) {
  reqinput <- c("y", "shape1", "shape2")
  checkNames1(input, reqinput)
  input <- input[reqinput]
  checkVector(input)
  
  for (i in seq_along(input)) {
    assign(names(input)[i], input[[i]])
  }
  if (any(!shape1>0)) stop("Parameter 'shape1' contains non-positive values.")
  if (any(is.infinite(shape1))) stop("Parameter 'shape1' contains infinite values.")
  if (any(!shape2>0)) stop("Parameter 'shape2' contains non-positive values.")
  if (any(is.infinite(shape2))) stop("Parameter 'shape2' contains infinite values.")
  
  return(list(y = y, shape1 = shape1, shape2 = shape2))
}
list_inputChecks$'beta' <- "check.beta"

################################################################################
### real line

### laplace
check.lapl <- function(input) {
  reqinput <- c("y", "location", "scale")
  checkNames1(input, reqinput)
  input <- input[reqinput]
  checkVector(input)
  
  for (i in seq_along(input)) {
    assign(names(input)[i], input[[i]])
  }
  if (any(is.infinite(location))) stop("Parameter 'location' contains infinite values.")
  if (any(!scale>0)) stop("Parameter 'scale' contains non-positive values.")
  if (any(is.infinite(scale))) stop("Parameter 'scale' contains infinite values.")
  
  return(list(y = y, location = location, scale = scale))
}
list_inputChecks$'laplace' <- "check.lapl"

### logistic
check.logis <- function(input) {
  reqinput <- c("y", "location", "scale")
  checkNames1(input, reqinput)
  input <- input[reqinput]
  checkVector(input)
  
  for (i in seq_along(input)) {
    assign(names(input)[i], input[[i]])
  }
  if (any(is.infinite(location))) stop("Parameter 'location' contains infinite values.")
  if (any(!scale>0)) stop("Parameter 'scale' contains non-positive values.")
  if (any(is.infinite(scale))) stop("Parameter 'scale' contains infinite values.")
  
  return(list(y = y, location = location, scale = scale))
}
list_inputChecks$'logistic' <- "check.logis"

### normal
check.norm <- function(input) {
  reqinput <- c("y", "mean", "sd")
  checkNames1(input, reqinput)
  input <- input[reqinput]
  checkVector(input)
  
  for (i in seq_along(input)) {
    assign(names(input)[i], input[[i]])
  }
  if (any(is.infinite(mean))) stop("Parameter 'mean' contains infinite values.")
  if (any(!sd>0)) stop("Parameter 'sd' contains non-positive values.")
  if (any(is.infinite(sd))) stop("Parameter 'sd' contains infinite values.")
  
  return(list(y = y, mean = mean, sd = sd))
}
list_inputChecks$'normal' <- "check.norm"

### normal-mixture
check.mixn <- function(input) {
  reqinput <- c("y", "m", "s", "w")
  checkNames1(input, reqinput)
  input <- input[reqinput]
  checkMatrix(input)
  
  for (i in seq_along(input)) {
    assign(names(input)[i], input[[i]])
  }
  if (any(is.infinite(m))) stop("Parameter 'm' contains infinite values.")
  if (any(!s>0)) stop("Parameter 's' contains non-positive values.")
  if (any(is.infinite(s))) stop("Parameter 's' contains infinite values.")
  if (any(w < 0 | w > 1)) stop("Parameter 'w' contains values not in [0, 1].")
  if (all.equal(apply(w, 1, sum), rep(1, dim(w)[1])) != TRUE) stop("Parameter 'w' contains weighting schemes which do not sum up to 1.")
  
  return(list(y = y, m = m, s = s, w = w))
}
list_inputChecks$'normal-mixture' <- "check.mixn"


### two-piece-normal
check.2pnorm <- function(input) {
  reqinput <- c("y", "m", "s1", "s2")
  checkNames1(input, reqinput)
  input <- input[reqinput]
  checkVector(input)
  
  for (i in seq_along(input)) {
    assign(names(input)[i], input[[i]])
  }
  if (any(is.infinite(m))) stop("Parameter 'm' contains infinite values.")
  if (any(!s1>0)) stop("Parameter 's1' contains non-positive values.")
  if (any(is.infinite(s1))) stop("Parameter 's1' contains infinite values.")
  if (any(!s2>0)) stop("Parameter 's2' contains non-positive values.")
  if (any(is.infinite(s2))) stop("Parameter 's2' contains infinite values.")
  
  return(list(y = y, m = m, s1 = s1, s2 = s2))
}
list_inputChecks$'two-piece-normal' <- "check.2pnorm"

### t
check.t <- function(input) {
  reqinput <- c("y", "df", "location", "scale")
  checkNames1(input, reqinput)
  input <- input[reqinput]
  checkVector(input)
  
  for (i in seq_along(input)) {
    assign(names(input)[i], input[[i]])
  }
  if (any(!df>0)) stop("Parameter 'df' contains non-positive values.")
  if (any(is.infinite(location))) stop("Parameter 'location' contains infinite values.")
  if (any(!scale>0)) stop("Parameter 'scale' contains non-positive values.")
  if (any(is.infinite(scale))) stop("Parameter 'scale' contains infinite values.")
  
  return(list(y = y, df = df, location = location, scale = scale))
}
list_inputChecks$'t' <- "check.t"


################################################################################
### non-negative

### exponential
check.exp <- function(input) {
  reqinput <- c("y", "rate")
  checkNames1(input, reqinput)
  input <- input[reqinput]
  checkVector(input)
  
  for (i in seq_along(input)) {
    assign(names(input)[i], input[[i]])
  }
  if (any(!rate>0)) stop("Parameter 'rate' contains non-positive values.")
  if (any(is.infinite(rate))) stop("Parameter 'rate' contains infinite values.")
  
  return(list(y = y, rate = rate))
}
list_inputChecks$'exponential' <- "check.exp"

### gamma
check.gamma <- function(input) {
  reqinput <- list(c("y", "shape", "rate"),
                   c("y", "shape", "scale")
  )
  choice <- checkNames2(input, reqinput)
  input <- input[reqinput[[choice]]]
  checkVector(input)
  
  for (i in seq_along(input)) {
    assign(names(input)[i], input[[i]])
  }
  if (any(!shape>0)) stop("Parameter 'shape' contains non-positive values.")
  if (any(!is.finite(shape))) stop("Parameter 'shape' contains infinite values.")
  if (choice == 1) {
    if (any(!rate>0)) stop("Parameter 'rate' contains non-positive values.")
    if (any(!is.finite(rate))) stop("Parameter 'rate' contains infinite values.")
    scale <- 1/rate
  } else if (choice == 2) {
    if (any(!scale>0)) stop("Parameter 'scale' contains non-positive values.")
    if (any(!is.finite(scale))) stop("Parameter 'scale' contains infinite values.")
  }
  
  return(list(y = y, shape = shape, scale = scale))
}
list_inputChecks$'gamma' <- "check.gamma"

### log-laplace
check.llapl <- function(input) {
  reqinput <- c("y", "locationlog", "scalelog")
  checkNames1(input, reqinput)
  input <- input[reqinput]
  checkVector(input)
  
  for (i in seq_along(input)) {
    assign(names(input)[i], input[[i]])
  }
  if (any(is.infinite(locationlog))) stop("Parameter 'locationlog' contains infinite values.")
  if (any(!scalelog>0)) stop("Parameter 'scalelog' contains non-positive values.")
  if (any(is.infinite(scalelog))) stop("Parameter 'scalelog' contains infinite values.")
  
  return(list(y = y, locationlog = locationlog, scalelog = scalelog))
}
list_inputChecks$'log-laplace' <- "check.llapl"

### log-logistic
check.llogis <- function(input) {
  reqinput <- c("y", "locationlog", "scalelog")
  checkNames1(input, reqinput)
  input <- input[reqinput]
  checkVector(input)
  
  for (i in seq_along(input)) {
    assign(names(input)[i], input[[i]])
  }
  if (any(is.infinite(locationlog))) stop("Parameter 'locationlog' contains infinite values.")
  if (any(!scalelog>0)) stop("Parameter 'scalelog' contains non-positive values.")
  if (any(is.infinite(scalelog))) stop("Parameter 'scalelog' contains infinite values.")
  
  return(list(y = y, locationlog = locationlog, scalelog = scalelog))
}
list_inputChecks$'log-logistic' <- "check.llogis"

### log-normal
check.lnorm <- function(input) {
  reqinput <- c("y", "meanlog", "sdlog")
  checkNames1(input, reqinput)
  input <- input[reqinput]
  checkVector(input)
  
  for (i in seq_along(input)) {
    assign(names(input)[i], input[[i]])
  }
  if (any(is.infinite(meanlog))) stop("Parameter 'meanlog' contains infinite values.")
  if (any(!sdlog>0)) stop("Parameter 'sdlog' contains non-positive values.")
  if (any(is.infinite(sdlog))) stop("Parameter 'sdlog' contains infinite values.")
  
  return(list(y = y, meanlog = meanlog, sdlog = sdlog))
}
list_inputChecks$'log-normal' <- "check.lnorm"

### truncated-normal
check.tn <- function(input) {
  reqinput <- c("y", "m", "s", "lb")
  checkNames1(input, reqinput)
  input <- input[reqinput]
  checkVector(input)
  
  for (i in seq_along(input)) {
    assign(names(input)[i], input[[i]])
  }
  if (any(is.infinite(m))) stop("Parameter 'm' contains infinite values.")
  if (any(!s>0)) stop("Parameter 's' contains non-positive values.")
  if (any(is.infinite(s))) stop("Parameter 's' contains infinite values.")
  if (any(is.infinite(lb))) stop("Parameter 'lb' contains infinite values.")
  
  return(list(y = y, m = m, s = s, lb = lb))
}
list_inputChecks$'truncated-normal' <- "check.tn"

################################################################################
### variable support

### gpd
check.gpd <- function(input) {
  reqinput <- c("y", "location", "scale", "shape")
  checkNames1(input, reqinput)
  input <- input[reqinput]
  checkVector(input)
  
  for (i in seq_along(input)) {
    assign(names(input)[i], input[[i]])
  }
  if (any(is.infinite(location))) stop("Parameter 'location' contains infinite values.")
  if (any(!scale>0)) stop("Parameter 'scale' contains non-positive values.")
  if (any(is.infinite(scale))) stop("Parameter 'scale' contains infinite values.")
  if (any(is.infinite(shape))) stop("Parameter 'shape' contains infinite values.")
  
  return(list(y = y, location = location, scale = scale, shape = shape))
}
list_inputChecks$'gpd' <- "check.gpd"

### gev
check.gev <- function(input) {
  reqinput <- c("y", "location", "scale", "shape")
  checkNames1(input, reqinput)
  input <- input[reqinput]
  checkVector(input)
  
  for (i in seq_along(input)) {
    assign(names(input)[i], input[[i]])
  }
  if (any(is.infinite(location))) stop("Parameter 'location' contains infinite values.")
  if (any(!scale>0)) stop("Parameter 'scale' contains non-positive values.")
  if (any(is.infinite(scale))) stop("Parameter 'scale' contains infinite values.")
  if (any(is.infinite(shape))) stop("Parameter 'shape' contains infinite values.")
      
  return(list(y = y, location = location, scale = scale, shape = shape))
}
list_inputChecks$'gev' <- "check.gev"