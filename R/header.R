# !!! Important !!! Set orientation = 1 for positive orientation, to - 1 for negative orientation
# (only used in functions that appear in the package: xx and xx_sample, where xx = (crps, logs, qs)
orientation <- -1

################################################################################
### xx_sample

crps_sample <- function(y, dat, method = "edf", w = NULL, bw = NULL, num_int = FALSE) {
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
        message("Used numerical integration for CPRS (tolerance = 1e-6).")
      }
    } else {
      stop("Unexpected choice of method - please select either edf or kde")
    }
    return(orientation * out)
  }

qs_sample <- function(y, dat, bw = NULL, num_int = FALSE) {
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
    message("Used numerical integration for qs computation (tolerance = 1e-6).")
    out <- 2*ff(y) - integrate(function(x) ff(x)^2, -Inf, Inf, rel.tol = 1e-6)$value
  }
  return(orientation * out)
}

logs_sample <- function(y, dat, bw = NULL) {
  input <- list(y = y, dat = dat)
  if (!is.null(bw)) {
    input$bw <- bw
  }
  check.sample(input)
  
  message(
    "Using the log score with kernel density estimation tends to be fragile -- see KLTG (2015) for details."
  )
  if (is.null(bw)) bw <- bw.nrd(dat)
  w <- rep(1 / length(dat), length(dat))
  s <- rep(bw, length(dat))
  out <- lsmixnC(w = w, m = dat, s = s, y = y)
  return(orientation * out)
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

crps <- function(y, family, ...) {
  ind <- c(match(family, names(list_crpsFunctions)),
           match(paste("crps.", family, sep=""), list_crpsFunctions)
  )
  ind <- ind[!is.na(ind)]
  ind <- unique(ind)
  ind <- names(list_crpsFunctions)[ind]
  if (length(ind) > 1) {
    stop("Ambiguous choice of parametric family - see details section of ?crps for a list of available choices.")
  } else if (length(ind) == 0) {
    stop("Could not find parametric family - see details section of ?crps for a list of available choices.")  
  }
  
  inputCheck <- eval(parse(text = list_inputChecks[ind]))
  input <- list(...)
  if (!missing(y)) input$y <- y
  input <- inputCheck(input)
  
  crps <- eval(parse(text = list_crpsFunctions[ind]))
  out <- do.call(crps, input)
  
  return(orientation*out)
}

qs <- function(y, family, ...) {
  ind <- c(match(family, names(list_qsFunctions)),
           match(paste("qs.", family, sep=""), list_qsFunctions)
  )
  ind <- ind[!is.na(ind)]
  ind <- unique(ind)
  ind <- names(list_qsFunctions)[ind]
  if (length(ind) > 1) {
    stop("Ambiguous choice of parametric family - see details section of ?qs for a list of available choices.")
  } else if (length(ind) == 0) {
    stop("Could not find parametric family - see details section of ?qs for a list of available choices.")  
  }
  
  inputCheck <- eval(parse(text = list_inputChecks[ind]))
  input <- list(...)
  if (!missing(y)) input$y <- y
  input <- inputCheck(input)
  
  qs <- eval(parse(text = list_qsFunctions[ind]))
  out <- do.call(qs, input)
  
  return(orientation*out)
}

logs <- function(y, family, ...) {
  ind <- c(match(family, names(list_lsFunctions)),
           match(paste("ls.", family, sep=""), list_lsFunctions)
  )
  ind <- ind[!is.na(ind)]
  ind <- unique(ind)
  ind <- names(list_lsFunctions)[ind]
  if (length(ind) > 1) {
    stop("Ambiguous choice of parametric family - see details section of ?logs for a list of available choices.")
  } else if (length(ind) == 0) {
    stop("Could not find parametric family - see details section of ?logs for a list of available choices.")  
  }
  
  inputCheck <- eval(parse(text = list_inputChecks[ind]))
  input <- list(...)
  if (!missing(y)) input$y <- y
  input <- inputCheck(input)
  
  logscore <- eval(parse(text = list_lsFunctions[ind]))
  out <- do.call(logscore, input)
  
  return(orientation*out)
}

################################################################################