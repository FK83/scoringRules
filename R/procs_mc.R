# Simulate Data from the Fox/West Model
sim.fw <- function(a, s, n, ndraws){
  # Compute parameters of marginal distribution (inverse Gamma)
  beta.ig <- c(0.5*(n+2), 0.5*n*s)
  # Compute parameters of corresponding t distribution
  beta.t <- c(beta.ig[2]/beta.ig[1], 2*beta.ig[1]) # variance, degrees of freedom
  # Simulate data
  sig2 <- rep(0, ndraws)
  sig2[1] <- 1
  psi <- 1/rgamma(ndraws, 0.5*(n+3), 0.5*n*s*(1-a^2))
  v <- a + sqrt( psi/(n*s) ) * rnorm(ndraws)
  for (j in 2:ndraws){
    sig2[j] <- psi[j] + (v[j]^2) * sig2[j-1]
  }  
  # Return list
  return(list(par.ig = beta.ig, sim = sig2))
}

# Make the ergodic dist for our simulation design
ergdist <- function(s, n){
  # Make t distribution parameters
  s <- sqrt(n*s/(n + 2))
  df <- n + 2
  # Compute properties of ergodic dist
  pdf <- function(z){
    sapply(z, function(l) dt(l/s, df)/s)
  }  
  cdf <- function(z){
    sapply(z, function(l) pt(l/s, df))
  } 
  list(pdf = pdf, cdf = cdf, s = s, df = df)
}

div.crps.t.mixn <- function(s.t, df.t, m.n, s.n) {
  cdf.t <- function(z) {
    pt(z/s.t, df.t)
  }
  cdf.n <- function(z) {
    sapply(z, function(s) mean(pnorm(s, mean = m.n, sd = s.n)))
  }
  div <- function(z) (cdf.t(z) - cdf.n(z))^2
  ret <- integrate(div, -Inf, Inf)
  return(ret$value)
}

div.crps.t.edf <- function(s.t, df.t, dat){
  dat <- sort(dat)
  F <- function(z) pt(z/s.t, df.t)
  m <- length(dat)
  out <- integrate(function(z) F(z)^2, -Inf, dat[1])$value
  for (j in 2:m){
    tmpf <- function(z) (F(z)-(j-1)/m)^2
    out <- out + integrate(tmpf, dat[j-1], dat[j])$value    
  } 
  out <- out + integrate(function(z) (1-F(z))^2, dat[m], Inf)$value  
  return(out) 
}

div.ls.t.mixn <- function(s.t, df.t, m.n, s.n, maxbd = 100){
  pdf.t <- function(z) {
    dt(z/s.t, df.t)/s.t
  }
  pdf.n <- function(z) {
    dmixnC(m = m.n, s = s.n, y = z)
  }
  f1 <- function(z) pdf.t(z)*(log(pdf.t(z)))
  f2 <- function(z) pdf.t(z)*(log(pdf.n(z)))
  
  # check boundaries of second function
  grid <- seq(from = -maxbd, to = maxbd, by = 1)
  check <- log(pdf.n(grid))
  lb <- grid[min(which(is.finite(check)))]
  ub <- grid[max(which(is.finite(check)))]
  
  ret1 <- integrate(f1, lb, ub)$value  
  ret2 <- integrate(f2, lb, ub)$value
  return(list(kl = ret1-ret2, lb = lb, ub = ub))
}

# convenience functions for log score
ls.n <- function(y, dat){
  m <- mean(dat)
  s <- sd(dat)
  mean(dnorm(y, mean = m, sd = s, log = TRUE))
}

make.bw <- function(dat, choice){
  if (choice == 1){
    s <- rep(bw.nrd(dat), length(dat))
  } else if (choice == 2){
    s <- rep(bw.ucv(dat), length(dat))
  } else if (choice == 3){
    s <- rep(bw.bcv(dat), length(dat))
  } else if (choice == 4){
    s <- rep(bw.SJ(dat), length(dat))
  }
  return(s)
}

# numerical integration for CRPS/mixture of normals
crps_mixn_int <- function(m, s, y, w = NULL){
  n <- length(m)
  if (is.null(w)){
    w <- rep(1/n, n)  
  } 
  Fmix <- function(z){
    sapply(z, function(r) sum(w*pnorm((r-m)/s)))
  }
  aux1 <- integrate(function(s) Fmix(s)^2, -Inf, y, rel.tol = 1e-6)$value
  aux2 <- integrate(function(s) (1 - Fmix(s))^2, y, Inf, rel.tol = 1e-6)$value
  (aux1 + aux2)
}

# scores
scores <- function(dat, m, s, y, sc = "crps", integration_threshold = 1000){
  
  # Sample size
  n_dat <- length(dat)
  
  # Some specifics, depending on wheter LS or CRPS is used
  if (sc == "crps"){
    out <- rep(0, 4)
    S1 <- crps
    S2 <- crps_sample
  } else if (sc == "logs"){
    out <- rep(0, 3)
    S1 <- logs
    S2 <- logs_sample
  }
  
  # Normal approx
  out[1] <- S1(y, family = "normal", mean = mean(dat), sd = sd(dat))
  
  # CKD
  if (sc == "crps" & n_dat > integration_threshold){
    out[2] <- crps_mixn_int(m = m, s = s, y = y)
  } else {
    out[2] <- S1(y, family = "normal-mixture", m = matrix(m, nrow = 1), s = matrix(s, nrow = 1),
                 w = matrix(rep(1/length(m), length(m)), nrow = 1))
  }
  
  # Kernel, ECDF
  if (sc == "crps"){
    if (n_dat < integration_threshold){
      out[3] <- S2(y, dat, method = "kde", show_messages = FALSE)  
    } else {
      out[3] <- S2(y, dat, method = "kde", num_int = TRUE, show_messages = FALSE)  
    }
    # Kernel
    
    # ECDF
    out[4] <- S2(y, dat, method = "edf", show_messages = FALSE)
  } else {
    # Kernel 
    out[3] <- S2(y, dat, show_messages = FALSE)
  }
  
  return(out)
}

n_helper <- function(x){
  ((x == 1000) + 2*(x == 5000) + 3*(x == 10000) + 4*(x == 20000) + 5*(x == 40000))
}


