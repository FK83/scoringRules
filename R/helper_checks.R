################################################################################
synonyms <- list(
  poisson = "pois",
  'negative-binomial' = "nbinom",
  uniform = "unif",
  laplace = "lapl",
  logistic = "logis",
  normal = "norm",
  'normal-mixture' = "mixnorm",
  'mixture-normal' = "mixnorm",
  'two-piece-exponential' = "2pexp",
  'two-piece-normal' = "2pnorm",
  exponential = "exp",
  'log-laplace' = "llapl",
  'log-logistic' = "llogis",
  'log-normal' = "lnorm"
)

# check existence of distribution family
#' @importFrom methods existsFunction
getFamily <- function(family, score) {
  family <- unique(family)
  ind <- match(family, names(synonyms), nomatch = 0L)
  family[ind > 0L] <- synonyms[ind]
  family <- unique(family)
  n <- length(family)
  
  if (n > 1) {
    stop(sprintf("Ambiguous choice of parametric family - see details section of ?%s for a list of available choices.",
                 paste0(score, ".numeric")))
  } else if (n == 0) {
    stop(sprintf("Could not find parametric family - see details section of ?%s for a list of available choices.",
                 paste0(score, ".numeric")))
  }
  if (!existsFunction(paste0(score, "_", family)) ||
      !existsFunction(paste0("check_", score, "_", family))) {
    stop(sprintf("Could not find parametric family - see details section of ?%s for a list of available choices.",
                 paste0(score, ".numeric")))
  }
  
  return(family)
}


################################################################################
### general parametric checks

# for only one choice of parameterization
checkNames1 <- function(required, given) {
  if (!all(required %in% given)) {
    stop(paste("Missing parameter.",
               paste("Given input:", paste(given, collapse=", ")),
               paste("Required input:", paste(required, collapse=", ")),
               sep="\n")
    )
  }
}

# for multiple parameterizations
checkNames2 <- function(required, given) {
  ind <- lapply(required, `%in%`, given)
  param.choice <- which(sapply(ind, all))
  if (length(param.choice) > 1) {
    stop(paste("Multiple parameterizations given. Please choose one.",
               paste("Given input:", paste(given, collapse = ", ")),
               paste("Required input:",
                     paste(lapply(required, paste, collapse = ", "),
                           collapse = "  OR  ")),
               sep="\n")
    ) 
  } else if (length(param.choice) == 0) {
    stop(paste("Missing parameter.",
               paste("Given input:", paste(given, collapse = ", ")),
               paste("Required input:",
                     paste(lapply(required, paste, collapse = ", "),
                           collapse = "  OR  ")),
               sep="\n")
    )
  }
}

checkNames3 <- function(optional, given) {
  ind <- match(optional, given, nomatch = 0)
  if (any(ind != 0)) {
    if (any(ind == 0)) {
      stop(
        paste(
          "Missing parameter.",
          paste("Given input:", paste(given, collapse = ", ")),
          paste("Optional input:", paste(optional, collapse = ", ")),
          "Using one of the optional parameters requires all of this group.",
          sep = "\n"
        )
      )
    } else {
      return(TRUE)
    }
  } else {
    return(FALSE)
  }
}

checkNumeric <- function(input, infinite_exception = NULL) {
  input_numeric <- sapply(input, is.numeric)
  if (any(!input_numeric)) {
    stop(paste("Non-numeric input:",
               paste(names(input)[!input_numeric], collapse = ", ")))
  }
  input_NA <- sapply(input, anyNA)
  if (any(input_NA)) {
    stop(paste("Input with missing values:",
               paste(names(input)[input_NA], collapse = ", ")))
  }
  input <- input[!names(input) %in% infinite_exception]
  input_infinite <- sapply(input, function(x) any(is.infinite(x)))
  if (any(input_infinite)) {
    stop(paste("Input with infinite values:",
               paste(names(input)[input_infinite], collapse = ", ")))
  }
}

# for scalar or vectorial parameters
checkVector <- function(input) {
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
  if (!is.vector(input$y)) stop("Non-scalar or non-vectorial input: y")
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

#### legacy checks ####

### logistic
check.logis <- function(input) {
  required <- c("y", "location", "scale")
  optional <- list(c("lower", "lmass"),
                   c("upper", "umass"))
  checkNames1(required, names(input))
  choice2 <- sapply(optional, checkNames3, names(input))
  
  used <- c(required, optional[choice2])
  input <- input[do.call(c, used)]
  
  checkVector(input)
  if (any(choice2)) {
    checkNumeric(input[!names(input) %in% c("lmass", "umass")])
    if (choice2[1] == TRUE) {
      if (is.character(input$lmass)){
        if (any(!input$lmass %in% c("trunc", "cens"))) {
          stop("Valid character input for 'lmass': 'trunc', 'cens'")
        }
      } else if (!is.numeric(input$lmass)) {
        stop("Type of 'lmass' is neither character nor numeric.")
      } else {
        if (any(input$lmass < 0 | input$lmass > 1)) {
          stop("Parameter 'lmass' contains values not in [0, 1].")
        }
      }
    }
    if (choice2[2] == TRUE) {
      if (is.character(input$umass)){
        if (any(!input$umass %in% c("trunc", "cens"))) {
          stop("Valid character input for 'lmass': 'trunc', 'cens'")
        }
      } else if (!is.numeric(input$umass)) {
        stop("Type of 'umass' is neither character nor numeric")
      } else {
        if (any(input$umass < 0 | input$umass > 1)) {
          stop("Parameter 'umass' contains values not in [0, 1].")
        }
      }
    }
    if (sum(choice2) == 2) {
      if (any(input$lower > input$upper)) {
        stop("Parameter 'lower' contains values greater than parameter 'upper'.")
      }  
    }
  } else {
    checkNumeric(input)
  }
  
  if (any(is.infinite(input$location))) stop("Parameter 'location' contains infinite values.")
  if (any(!input$scale>0)) stop("Parameter 'scale' contains non-positive values.")
  if (any(is.infinite(input$scale))) stop("Parameter 'scale' contains infinite values.")

  return(input)
}

### normal
check.norm <- function(input) {
  required <- list(c("y", "mean", "sd"),
                   c("y", "location", "scale"))
  optional <- list(c("lower", "lmass"),
                   c("upper", "umass"))
  choice <- checkNames2(required, names(input)) # object type: integer
  choice2 <- sapply(optional, checkNames3, names(input)) # object type: logical
  
  used <- c(required[choice], optional[choice2]) # object type: list
  input <- input[do.call(c, used)]
  
  checkVector(input)
  if (any(choice2)) {
    checkNumeric(input[!names(input) %in% c("lmass", "umass")])
    if (choice2[1] == TRUE) {
      if (is.character(input$lmass)){
        if (any(!input$lmass %in% c("trunc", "cens"))) {
          stop("Valid character input for 'lmass': 'trunc', 'cens'")
        }
      } else if (!is.numeric(input$lmass)) {
        stop("Type of 'lmass' is neither character nor numeric.")
      } else {
        if (any(input$lmass < 0 | input$lmass > 1)) {
          stop("Parameter 'lmass' contains values not in [0, 1].")
        }
      }
    }
    if (choice2[2] == TRUE) {
      if (is.character(input$umass)){
        if (any(!input$umass %in% c("trunc", "cens"))) {
          stop("Valid character input for 'lmass': 'trunc', 'cens'")
        }
      } else if (!is.numeric(input$umass)) {
        stop("Type of 'umass' is neither character nor numeric")
      } else {
        if (any(input$umass < 0 | input$umass > 1)) {
          stop("Parameter 'umass' contains values not in [0, 1].")
        }
      }
    }
    if (sum(choice2) == 2) {
      if (any(input$lower > input$upper)) {
        stop("Parameter 'lower' contains values greater than parameter 'upper'.")
      }  
    }
  } else {
    checkNumeric(input)
  }
  
  if (choice == 1) {
    if (any(is.infinite(input$mean))) stop("Parameter 'mean' contains infinite values.")
    if (any(!input$sd>0)) stop("Parameter 'sd' contains non-positive values.")
    if (any(is.infinite(input$sd))) stop("Parameter 'sd' contains infinite values.")
    names(input)[names(input) == "mean"] <- "location"
    names(input)[names(input) == "sd"] <- "scale"
  } else if (choice == 2) {
    if (any(is.infinite(input$location))) stop("Parameter 'location' contains infinite values.")
    if (any(!input$scale>0)) stop("Parameter 'scale' contains non-positive values.")
    if (any(is.infinite(input$scale))) stop("Parameter 'scale' contains infinite values.")
  }
  
  return(input)
}

### t
check.t <- function(input) {
  required <- c("y", "df", "location", "scale")
  optional <- list(c("lower", "lmass"),
                   c("upper", "umass"))
  checkNames1(required, names(input))
  choice2 <- sapply(optional, checkNames3, names(input))
  
  used <- c(required, optional[choice2])
  input <- input[do.call(c, used)]
  
  checkVector(input)
  if (any(choice2)) {
    checkNumeric(input[!names(input) %in% c("lmass", "umass")])
    if (choice2[1] == TRUE) {
      if (is.character(input$lmass)){
        if (any(!input$lmass %in% c("trunc", "cens"))) {
          stop("Valid character input for 'lmass': 'trunc', 'cens'")
        }
      } else if (!is.numeric(input$lmass)) {
        stop("Type of 'lmass' is neither character nor numeric.")
      } else {
        if (any(input$lmass < 0 | input$lmass > 1)) {
          stop("Parameter 'lmass' contains values not in [0, 1].")
        }
      }
    }
    if (choice2[2] == TRUE) {
      if (is.character(input$umass)){
        if (any(!input$umass %in% c("trunc", "cens"))) {
          stop("Valid character input for 'lmass': 'trunc', 'cens'")
        }
      } else if (!is.numeric(input$umass)) {
        stop("Type of 'umass' is neither character nor numeric")
      } else {
        if (any(input$umass < 0 | input$umass > 1)) {
          stop("Parameter 'umass' contains values not in [0, 1].")
        }
      }
    }
    if (sum(choice2) == 2) {
      if (any(input$lower > input$upper)) {
        stop("Parameter 'lower' contains values greater than parameter 'upper'.")
      }  
    }
  } else {
    checkNumeric(input)
  }
  
  if (any(!input$df>0)) stop("Parameter 'df' contains non-positive values.")
  if (any(is.infinite(input$location))) stop("Parameter 'location' contains infinite values.")
  if (any(!input$scale>0)) stop("Parameter 'scale' contains non-positive values.")
  if (any(is.infinite(input$scale))) stop("Parameter 'scale' contains infinite values.")
  
  return(input)
}

### normal-mixture
check.mixnorm <- function(input) {
  required <- c("y", "m", "s", "w")
  optional <- c("exact", "rel.tol")
  checkNames1(required, names(input))
  checkMatrix(input[required])
  ind <- sapply(optional, checkNames3, names(input))
  input <- input[c(required, optional[ind])]
  if (ind[2] == TRUE) {
    if (length(input$rel.tol) != 1) {
      stop("Parameter 'rel.tol' needs to be of length 1.")
    }
  }
  if (ind[1] == TRUE) {
    if (length(input$exact) != 1) {
      stop("Parameter 'exact' needs to be of length 1.")
    }
    if (!input$exact %in% c(TRUE, FALSE)) {
      stop("Parameter 'exact' needs to be TRUE or FALSE.")
    }
    checkNumeric(input[names(input) != "exact"])
  } else {
    checkNumeric(input)
  }
  
  if (any(is.infinite(input$m))) stop("Parameter 'm' contains infinite values.")
  if (any(!input$s>0)) stop("Parameter 's' contains non-positive values.")
  if (any(is.infinite(input$s))) stop("Parameter 's' contains infinite values.")
  if (any(input$w < 0 | input$w > 1)) stop("Parameter 'w' contains values not in [0, 1].")
  if (all.equal(apply(input$w, 1, sum), rep(1, dim(input$w)[1])) != TRUE) stop("Parameter 'w' contains weighting schemes which do not sum up to 1.")
  
  return(input)
}
