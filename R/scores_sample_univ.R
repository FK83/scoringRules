#' Scoring Rules for Simulated Forecast Distributions
#' 
#' Calculate scores (CRPS, LogS, DSS) given observations and draws from the predictive distributions.
#' 
#' @param y vector of realized values.
#' @param dat vector or matrix (depending on \code{y}; see details)
#'  of simulation draws from forecast distribution. 
#' @param method string; approximation method. Options:
#'  "edf" (empirical distribution function) and
#'  "kde" (kernel density estimation).
#' @param w optional; vector or matrix (matching \code{dat}) of weights for method \code{"edf"}.
#' @param bw optional; vector (matching \code{y}) of bandwidths for kernel density
#' estimation; see details.
#' @param num_int logical; if TRUE numerical integration is used for method \code{"kde"}.
#' @param show_messages logical; display of messages (does not affect
#'  warnings and errors).
#'  
#' @return
#' Value of the score. \emph{A lower score indicates a better forecast.}
#' 
#' @references
#' \emph{Evaluating simulation based forecast distributions:}
#' 
#' Krueger, F., Lerch, S., Thorarinsdottir, T.L. and T. Gneiting (2020): `Predictive inference based on Markov chain Monte Carlo output', \emph{International Statistical Review}, forthcoming. \doi{10.1111/insr.12405}
#'  
#' \emph{Empirical quantile decomposition of the CRPS:}
#'  
#' Laio, F. and S. Tamea (2007):
#' `Verification tools for probabilistic forecasts of continuous
#' hydrological variables',
#' Hydrology and Earth System Sciences, 11, 1267-1277. \doi{10.5194/hess-11-1267-2007}
#'  
#' @author Alexander Jordan, Fabian Krueger, Sebastian Lerch
#'  
#' @details 
#' For a vector \code{y} of length n, \code{dat} should be given as a matrix
#' with n rows. If \code{y} has length 1, then \code{dat} may be a vector.
#' 
#' \code{\link{crps_sample}} employs an empirical version of the quantile
#' decomposition of the CRPS (Laio and Tamea, 2007) when using
#' \code{method = "edf"}. For \code{method = "kde"}, it uses kernel density
#' estimation using a Gaussian kernel. The logarithmic score always uses kernel density estimation.
#' 
#' The bandwidth (\code{bw}) for kernel density estimation can be
#' specified manually, in which case it must be a positive number. If
#' \code{bw == NULL}, the bandwidth is selected using the core function
#' \code{\link{bw.nrd}}. Numerical integration may speed up computation for
#' \code{\link{crps_sample}} in case of large samples \code{dat}.
#'  
#' @examples
#' \dontrun{
#' 
#' # y has length greater than 1
#' y <- 1:2
#' sample <- matrix(rnorm(20), nrow = 2)
#' crps_sample(y = y, dat = sample)
#' logs_sample(y = y, dat = sample)
#' 
#' y <- 1:2
#' sample <- rnorm(10)
#' crps_sample(y = y, dat = sample) # error
#' 
#' # y has length 1
#' y <- 1
#' sample <- rnorm(10)
#' crps_sample(y = y, dat = sample)
#' 
#' sample <- matrix(rnorm(10), nrow = 1)
#' crps_sample(y = y, dat = sample)
#' 
#' sample <- matrix(rnorm(20), nrow = 2)
#' crps_sample(y = y, dat = sample) # error
#' }
#' 
#' @name scores_sample_univ
#' @importFrom stats bw.nrd
NULL

#' @rdname scores_sample_univ
#' @export
crps_sample <- function(y, dat, method = "edf", w = NULL, bw = NULL, 
                        num_int = FALSE, show_messages = TRUE) {
  input <- list(y = y, dat = dat)
  if (method == "edf") {
    if (show_messages) {
      if (!is.null(bw)) message("Parameter 'bw' is ignored for edf method.")
      if (num_int) message("Parameter 'num_int' is ignored for edf method.")
    }
    if (!is.null(w)) input$w <- w
    if (identical(length(y), 1L) && is.vector(dat)) {
      check_sample(input)
      crps_edf(y, dat, w)
    } else {
      check_sample2(input)
      sapply(seq_along(y),
             function(i) crps_edf(y[i], dat[i, ], w[i, ]))
    }
  } else if (method == "kde") {
    if (show_messages) {
      if (!is.null(w)) message("Parameter 'w' is ignored for kde method.")
      if (num_int) message("Used numerical integration for CRPS computation (tolerance = 1e-6).")
    }
    if (!is.null(bw)) input$bw <- bw
    if (identical(length(y), 1L) && is.vector(dat)) {
      check_sample(input)
      crps_kdens(y, dat, bw, num_int)
    } else {
      check_sample2(input)
      crps_kdens(y, dat, bw, num_int)
    }
  } else {
    stop("Unexpected choice of method - please select either 'edf' or 'kde'.")
  }
}


#' @rdname scores_sample_univ
#' @export
logs_sample <- function(y, dat, bw = NULL, show_messages = FALSE) {
  input <- list(y = y, dat = dat)
  input$bw <- bw
  
  if (show_messages)
    message("Using the log score with kernel density estimation tends to be fragile -- see KLTG (2019) for details.")
  if (identical(length(y), 1L) && is.vector(dat)) {
    check_sample(input)
    if (is.null(bw)) bw <- bw.nrd(dat)
    n <- length(dat)
    w <- rep(1 / n, n)
    s <- rep(bw, n)
    lsmixnC(w, dat, s, y)
  } else {
    check_sample2(input)
    if (is.null(bw)) bw <- apply(dat, 1, bw.nrd)
    w <- rep(1, ncol(dat))
    s <- matrix(bw, nrow(dat), ncol(dat))
    sapply(seq_along(y), function(i) lsmixnC(w, dat[i, ], s[i, ], y[i]))
  }
}


#' @rdname scores_sample_univ
#' @export
dss_sample <- function(y, dat, w = NULL) {
  input <- list(y = y, dat = dat)
  if (!is.null(w)) input$w <- w
  if (identical(length(y), 1L) && is.vector(dat)) {
    check_sample(input)
    dss_edf(y, dat, w)
  } else {
    check_sample2(input)
    sapply(seq_along(y),
           function(i) dss_edf(y[i], dat[i, ], w[i, ]))
  }
}

#### helper functions ####

# (weighted) empirical distribution
crps_edf <- function(y, dat, w = NULL) {
  if (is.null(w)) {
    c_1n <- 1 / length(dat)
    x <- sort(dat)
    a <- seq.int(0.5 * c_1n, 1 - 0.5 * c_1n, length.out = length(dat))
    f <- function(s) 2 * c_1n * sum(((s < x) - a) * (x - s))
  } else {
    if (!identical(length(dat), length(w)) || any(w < 0, na.rm = TRUE)) {
      return(rep(NaN, length(y)))
    }
    ord <- order(dat)
    x <- dat[ord]
    w <- w[ord]
    p <- cumsum(w)
    P <- p[length(p)]
    a <- (p - 0.5 * w) / P
    f <- function(s) 2 / P * sum(w * ((s < x) - a) * (x - s))
  }
  sapply(y, f)
}

dss_edf <- function(y, dat, w = NULL) {
  if (is.null(w)) {
    m <- mean(dat)
    v <- mean(dat^2) - m^2
  } else {
    if (!identical(length(dat), length(w)) || any(w < 0, na.rm = TRUE)) {
      return(rep(NaN, length(y)))
    }
    W <- sum(W)
    m <- sum(w * dat) / W
    v <- sum(w * dat^2) / W - m^2
  }
  sapply(y, function(s) (s - m)^2 / v + log(v))
}

# kernel density estimation
crps_kdens <- function(y, dat, bw = NULL, num_int = FALSE) {
  if (is.vector(dat)) {
    n <- length(dat)
    if (is.null(bw)) bw <- bw.nrd(dat)
    dim(dat) <- c(1L, n)
    s <- matrix(bw,  nrow = 1L, ncol = n)
    w <- matrix(1/n, nrow = 1L, ncol = n)
    if (num_int == FALSE) {
      sapply(y, function(x) crps_mixnorm(x, dat, s, w))
    } else {
      sapply(y, function(x) crps_mixnorm_int(x, dat, s, w))
    }
  } else if (is.matrix(dat)) {
    n1 <- dim(dat)[1L]
    n2 <- dim(dat)[2L]
    if (is.null(bw)) bw <- apply(dat, 1L, bw.nrd)
    s <- matrix(bw,   nrow = n1, ncol = n2)
    w <- matrix(1/n2, nrow = n1, ncol = n2)
    if (num_int == FALSE) {
      crps_mixnorm(y, dat, s, w)
    } else {
      crps_mixnorm_int(y, dat, s, w)
    }
  } else {
    stop("Invalid data type for 'dat'.")
  }
}


#### input checks for sample functions ####
check_sample <- function(input) {
  checkNumeric(input)
  
  input_isvector <- sapply(input, is.vector)
  if (!all(input_isvector)) {
    stop(paste("Non-vector input:",
               paste(names(input)[!input_isvector], collapse=", ")))
  }
  
  input_lengths <- sapply(input, length)
  max_length <- max(input_lengths)
  ref_lengths <- c(y = 1L, dat = max_length, w = max_length, bw = 1L)
  if (!identical(input_lengths, ref_lengths[names(input)])) {
    ref_lengths2 <- c(y = "1", dat = "n", w = "n", bw = "1")
    stop(
      paste(
        "Incompatible input vector lengths.",
        sprintf("Lengths of (%s) should be (%s).",
                paste(names(input), collapse = ", "),
                paste(ref_lengths2[names(input)], collapse = ", ")),
        sprintf("Given lengths: %s", paste(input_lengths, collapse = ", ")),
        sep = "\n")
    )
  }
  
  if (!is.null(input$w)) {
    if (any(input$w < 0)) {
      stop("Weight parameter 'w' contains negative values.")
    }
  }
  if (!is.null(input$bw)) {
    if (input$bw < 0) {
      stop("Bandwidth parameter 'bw' is negative.")
    }
  }
}

check_sample2 <- function(input) {
  checkNumeric(input)
  
  input_isvector <- sapply(input[names(input) %in% c("y", "bw")], is.vector)
  input_ismatrix <- sapply(input[names(input) %in% c("dat", "w")], is.matrix)
  input_dim1 <- sapply(input, function(x) {
    if (is.vector(x)) length(x) else dim(x)[1L]
  })
  input_dim2 <- sapply(input[names(input) %in% c("dat", "w")], function(x) {
    if (is.matrix(x)) dim(x)[2L]
  })
  
  if (!all(input_isvector) ||
      !all(input_ismatrix) ||
      !identical(length(unique(input_dim1)), 1L) ||
      !identical(length(unique(input_dim2)), 1L)) {
    
    reference_formats <- c(
      y = "vector[1:n]",
      dat = "matrix[1:n, 1:m]",
      w = "matrix[1:n, 1:m]",
      bw = "vector[1:n]"
    )
    input_formats <- sapply(input, function(x) {
      if (is.vector(x)) {
        sprintf("vector[1:%i]", length(x))
      } else if (is.matrix(x)) {
        sprintf("matrix[1:%i, 1:%i]", dim(x)[1L], dim(x)[2L])
      } else {
        "unidentified"
      }
    })
    stop(
      paste(
        "Incompatible input objects.",
        sprintf("Expected input for (%s): %s.",
                paste(names(input), collapse = ", "),
                paste(reference_formats[names(input)], collapse = ", ")),
        sprintf("   Given input for (%s): %s.",
                paste(names(input), collapse = ", "),
                paste(input_formats, collapse = ", ")),
        sep = "\n"
      )
    )
  }
  
  if (!is.null(input$w)) {
    if (any(input$w < 0)) {
      stop("Weight parameter 'w' contains negative values.")
    }
  }
  if (!is.null(input$bw)) {
    if (input$bw < 0) {
      stop("Bandwidth parameter 'bw' is negative.")
    }
  }
}
