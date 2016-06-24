# !!! Important !!! Set orientation = 1 for positive orientation, to - 1 for negative orientation
# (only used in functions that appear in the package: xx and xx_sample, where xx = (crps, logs, qs)
orientation <- -1

#' @export crps logs crps_sample logs_sample flapl f2pexp fmixnorm f2pnorm ft fllapl fllogis 
#' @importFrom Rcpp evalCpp
#' @useDynLib scoringRules

################################################################################
### xx_sample

crps_sample <- function(y, dat, method = "edf", w = NULL, bw = NULL, 
                        num_int = FALSE, show_messages = TRUE) {
    input <- list(y = y, dat = dat)
    if (!is.null(w)) {
      input$w <- w
    }
    if (!is.null(bw)) {
      input$bw <- bw
    }
    check.sample(input)
    
    # Throw error if weights are not equal
    if (method == "kde" & !is.null(w)) {
      warning("Observation weights not implemented for kernel density estimation - edf used instead.")
      method <- "edf"
    }
    
    if (method == "edf") {
      out <- crps.edf(dat = dat, y = y, w = w)
    } else if (method == "kde") {
      # Bandwidth
      if (is.null(bw)) {
        bw <- bw.nrd(dat)
      }
      if (num_int == FALSE) {
        # Exact formula
        out <- crps.kdens(dat = dat, y = y, bw = bw)
      } else {
        # Numerical integration
        FF <- function(x) {
          sapply(x, function(zz) mean(pnorm(zz, mean = dat, sd = bw)))
        }
        aux1 <- integrate(function(s) FF(s)^2,-Inf, y, rel.tol = 1e-6)$value
        aux2 <- integrate(function(s) (1 - FF(s))^2, y, Inf, rel.tol = 1e-6)$value
        out <- -(aux1 + aux2)
        # Message
        if (show_messages) 
          message("Used numerical integration for CRPS computation (tolerance = 1e-6).")
      }
    } else {
      stop("Unexpected choice of method - please select either edf or kde")
    }
    return(out)
  }

qs_sample <- function(y, dat, bw = NULL, num_int = FALSE, show_messages = TRUE) {
  input <- list(y = y, dat = dat)
  if (!is.null(bw)) {
    input$bw <- bw
  }
  check.sample(input)
  
  if (is.null(bw)) {
    bw <- bw.nrd(dat)
  }
  if (num_int == FALSE) {
    w <- rep(1 / length(dat), length(dat))
    s <- rep(bw, length(dat))
    out <- qsmixnC(w = w, m = dat, s = s, y = y)
  } else {
    ff <- function(x) {
      sapply(x, function(zz) mean(dnorm(zz, mean = dat, sd = bw)))
    }
    if (show_messages) 
      message("Used numerical integration for qs computation (tolerance = 1e-6).")
    out <- 2*ff(y) - integrate(function(x) ff(x)^2, -Inf, Inf, rel.tol = 1e-6)$value
  }
  return(out)
}

logs_sample <- function(y, dat, bw = NULL, show_messages = TRUE) {
  input <- list(y = y, dat = dat)
  if (!is.null(bw)) {
    input$bw <- bw
  }
  check.sample(input)
  
  if (show_messages)
    message("Using the log score with kernel density estimation tends to be fragile -- see KLTG (2016) for details.")
  if (is.null(bw)) bw <- bw.nrd(dat)
  w <- rep(1 / length(dat), length(dat))
  s <- rep(bw, length(dat))
  out <- lsmixnC(w = w, m = dat, s = s, y = y)
  return(out)
}

# input checks for sample functions
check.sample <- function(input) {
  
  input_isnumeric <- sapply(input, is.numeric)
  if (!all(input_isnumeric)) {
    stop(paste("Non-numeric input:", paste(names(input)[!input_isnumeric], collapse=", ")))
  }
  
  input_isvector <- sapply(input, is.vector)
  if (!all(input_isvector)) {
    stop(paste("Non-scalar or non-vectorial input:", paste(names(input)[!input_isvector], collapse=", ")))
  }
  
  input_lengths <- sapply(input, length)
  max_length <- max(input_lengths)
  ref_length <- do.call(c, list(y = 1, dat = max_length, w = max_length, bw = 1))
  ref_length2 <- do.call(c, list(y = "1", dat = "n", w = "n", bw = "1"))
  length_diffs <- input_lengths - ref_length[names(input)]
  if (any(length_diffs != 0)) {
    stop(paste("Incompatible input vector lengths.",
               sprintf("Lengths of (%s) should be (%s).",
                       paste(names(input), collapse = ", "),
                       paste(ref_length2[names(input)], collapse = ", ")
               ),
               sprintf("Given lengths: %s", paste(input_lengths, collapse = ", ")),
               sep = "\n")
    )
  }
  
  if (!is.null(input$w)) {
    w <- input$w
    if (any(w < 0 | w > 1)) {
      stop("Weight parameter 'w' contains values not in [0, 1].")
    }
    if (!isTRUE(all.equal(sum(w), 1))) {
      stop("Weight parameter 'w' does not sum up to 1.")
    }
  }
  if (!is.null(input$bw)) {
    if (input$bw < 0) {
      stop("Bandwidth parameter 'bw' is negative.")
    }
  }
}


################################################################################
### parametric

checkFamily <- function(family, score) {
  family <- unique(family)
  ind <- match(family, names(synonyms), nomatch = 0)
  family[ind > 0] <- synonyms[ind]
  family <- unique(family)
  n <- length(family)
  
  if (n > 1) {
    stop(sprintf("Ambiguous choice of parametric family - see details section of ?%s for a list of available choices.",
                 score))
  } else if (n == 0) {
    stop(sprintf("Could not find parametric family - see details section of ?%s for a list of available choices.",
                 score))
  }
  if (!existsFunction(paste0(score, ".", family)) | !existsFunction(paste0("check.", family))) {
    stop(sprintf("Could not find parametric family - see details section of ?%s for a list of available choices.",
                 score))
  }
  
  return(family)
}

crps <- function(y, family, ...) {
  family <- checkFamily(family, "crps")
  checkInput <- get(paste0("check.", family))
  calculateCRPS <- get(paste0("crps.", family))
  
  if (is.list(y)) {
    input <- c(y, list(...))
  } else {
    input <- list(y = y, ...)
  }
  
  input <- checkInput(input)
  out <- do.call(calculateCRPS, input)
  
  if (any(out < 0)) {
    warning("Negative CRPS values. Check parameter combinations and/or contact package maintainer(s).")
  }
  if (any(is.na(out))) {
    warning("Missing CRPS values. Probably due to numerical instabilities as a result of extreme parameter choices.")
  }
  
  return(out)
}

qs <- function(y, family, ...) {
  family <- checkFamily(family, "qs")
  checkInput <- get(paste0("check.", family))
  calculateQS <- get(paste0("qs.", family))
  
  if (is.list(y)) {
    input <- c(y, list(...))
  } else {
    input <- list(y = y, ...)
  }
  
  input <- checkInput(input)
  out <- do.call(calculateQS, input)
  
  return(out)
}

logs <- function(y, family, ...) {
  family <- checkFamily(family, "ls")
  checkInput <- get(paste0("check.", family))
  calculateLS <- get(paste0("ls.", family))
  
  if (is.list(y)) {
    input <- c(y, list(...))
  } else {
    input <- list(y = y, ...)
  }
  
  input <- checkInput(input)
  out <- do.call(calculateLS, input)
  
  return(out)
}

################################################################################