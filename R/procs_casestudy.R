#' @importFrom Rcpp evalCpp
#' @importFrom MASS mvrnorm
#' @importFrom graphics axis legend lines mtext plot points segments title
#' @importFrom stats aggregate bw.SJ bw.nrd bw.ucv bw.bcv dnorm dt integrate pnorm pt quantile rchisq rgamma rnorm runif sd

############################################################################################################################
# Code for estimating (and forecasting with) an autoregressive model with Markov Switching parameters
# Adapted from Matlab code supplied by Gianni Amisano via his website at University of Brescia (web page now defunct)
############################################################################################################################

# Function to draw from conditional posterior distribution of beta
drawBeta <- function(yv, x1, x2, hv, sv, T, k1, k2, nm, betaPriorPrec, betaPriorPrecMean){
  nreg <- k1 + k2*nm
  scalev <- sqrt(hv[sv])
  if (k2 > 0){
    w2 <- matrix(sv, length(sv), k2) == 1
    for (im in 2:nm){
      w2 <- cbind(w2, matrix(sv, length(sv), k2) == im)
    }
    w2 <- matrix(x2, nrow(x2), nm*ncol(x2)) * w2 
  } else {
    w2 <- {}
  }
  x1w2 <- cbind(x1, w2)
  w <- matrix(scalev, length(scalev), nreg) * x1w2
  wv <- scalev * yv
  betaPostVar <- solve(t(w)%*% w + betaPriorPrec)
  betaPostMean <- betaPostVar %*% (t(w) %*% wv + betaPriorPrecMean)
  betaDraw <- mvrnorm(mu = betaPostMean, Sigma = betaPostVar)
  epsv <- yv - x1w2 %*% betaDraw
  return(list(betaDraw = betaDraw, epsv = epsv))
}

# Simplified version of function "draw_dirichlet()"
# M = sample size, r = Dirichlet parameter (length = dim of random variable)
drawDirichlet <- function(M, r){
  y <- t(matrix(rchisq(M*length(r), rep(r, M)), length(r), M))
  return(y/rowSums(y))
}

# Analogue to matlab function MSLinRegDrawP.m
drawP <- function(T, nm, sv, rm){
	Pdraw <- matrix(0, nm, nm)
	sm <- cbind(sv[1:(T-1)], sv[2:T])
	tm <- matrix(0, nm, nm)
	for (im in 1:nm){
		for (jm in 1:nm){
			tm[im, jm] <- sum( (sm[, 1] == im) & (sm[, 2] == jm) )                
		}
		Pdraw[im, ] <- drawDirichlet(1, rm[im, ] + tm[im, ])   
	}
	pv <- rbind(t(Pdraw) - diag(nm), matrix(1, 1, nm))
	test <- det(t(pv) %*% pv)
	if (test != 0){ 
		pv <- solve(t(pv) %*% pv) %*% t(pv);
		pv <- matrix(pv[, nm + 1])           
	} else {
		pv <- matrix(rep(1/nm, nm))
	}
	return(list(Pdraw = Pdraw, pv = pv, tm = tm))
}

# Draw from multinomial dist
drawMultinom <- function(ndraw, prob){
  prob <- cumsum(prob)
  aux <- sapply(runif(ndraw), function(s) s < prob)
  return(apply(aux, 2, function(s) which.max(s)))
}

# Draw precisions
drawHv <- function(T, nm, indSigma, Sv, nuv, epsv, sv){
  hvDraw <- rep(0, nm)
  if (indSigma != 0){
    for (im in 1:nm){
      epiv <- epsv[sv == im]
      S1 <- Sv[im] + sum(epiv^2)
      hvDraw[im] <- rchisq(1, nuv[im] + length(epiv)) / S1      
    }    
  } else {
    S1 <- Sv[1] + sum(epsv^2)
    hvDraw <- rep(1, nm) * rchisq(1, nuv[1] + T) / S1    
  }
  return(hvDraw)  
}

# filtering -> maybe write this in C++
filterMarkovMixture <- function(p, P, lnpdat, indSS){
  T <- nrow(lnpdat)
  m <- ncol(lnpdat)
  simstate <- 0
  filprob <- matrix(0, T, m)
  lnl <- rep(0, T)
  pit1t1 <- p
  for (t in 1:T){
    pitt1 <- t(P) %*% pit1t1  
    lnpmax <- max(lnpdat[t, pitt1 > 0])
    pitt <- pitt1 * matrix(exp(lnpdat[t, ] - lnpmax))
    lkl <- sum(pitt)
    if (lkl == 0){
      print(pitt1)
    }
    pitt <- pitt / lkl
    lnl[t] <- log(lkl) + lnpmax
    filprob[t, ] <- t(pitt)
    pit1t1 <- pitt    
  }
  if (indSS != 0){
    simstate <- rep(0, T)
    simstate[T] <- drawMultinom(1, t(filprob[T, ]))
    for (t in (T-1):1){
      p1 <- t(filprob[t, ]) * P[, simstate[t+1]]
      p1 <- p1 / sum(p1)
      simstate[t] <- drawMultinom(1, p1)      
    }   
  }
 return(list(lnl = lnl, filprob = filprob, simstate = simstate))   
}

# helper function to multiply a matrix by itself, mult times
matmult <- function(mat, mult){
  if (mult == 0){
    out <- diag(nrow(mat))
  } else if (mult == 1){
    out <- mat
  } else {
    out <- mat
    for (c in 1:(mult-1)) out <- out %*% mat    
  }
  return(out)
}

# draw state variable
drawState <- function(T, k1, k2, nm, indSigma, p, P, yv, x1, x2, beta, hv, useC = TRUE){
  b1 <- beta[1:k1]
  b2 <- beta[(k1+1) : (k1 + k2*nm)]
  meanmat <- rep(0, T, nm)
  if (k1 != 0){
    tmp <- x1 %*% b1
    meanmat <- meanmat + matrix(tmp, nrow(tmp), nm * ncol(tmp))
  }
  if (k2 != 0){
    b2 <- matrix(b2, k2, nm)
    meanmat <- meanmat + x2 %*% b2
  }
  sigm <- t(matrix(1/sqrt(hv), length(hv), T))  
  aux <- matrix(yv, length(yv), nm)
  lnpdat <- log(sapply(1:nm, function(z) dnorm(aux[, z], mean = meanmat[, z], sd = sigm[, z])))
  if (useC == TRUE){
    aux <- filterMarkovMixtureC(p, P, lnpdat)
  } else {
    aux <- filterMarkovMixture(p, P, lnpdat, 1)
  }
  return(list(lnl = aux$lnl, filprob = aux$filprob, simstate = aux$simstate))  
}

# predictive density
predDensAR <- function(y, p, theta, nm, fpv, beta.switch, variance.switch, hmax = 10){
  
  # Extract coefficients
  if (beta.switch == TRUE){
    indb <- nm*(p+1) 
    betas <- matrix(0, p+1, nm)
    for (l in 1:nm) betas[, l] <- theta[((l-1)*(p+1)+1):(l*(p+1))]
  } else {
    indb <- p + 1
    betas <- matrix(rep(theta[1:(p+1)], nm), p+1, nm)   
  }  
  if (variance.switch == TRUE){
    indv <- (indb+nm)
    vs <- 1/theta[(indb+1):(indv)]
  } else { 
    indv <- indb + 1
    vs <- rep(1/theta[indv], nm)
  } 
  P <- matrix(theta[(indv + 1): length(theta)], nm, nm)
  
  # Simulate future states
  fs <- c(drawMultinomC(t(P) %*% fpv), rep(0, hmax - 1))
  for (hh in 2:hmax){
    aux <- rep(0, nm)
    aux[fs[hh-1]] <- 1
    fs[hh] <- drawMultinomC(t(P) %*% aux)    
  }
  
  # helper function
  replace.reg <- function(x, fc){
    if (length(x) == 1){
      return(fc)
    } else {
      return(cbind(fc, x[-length(x)]))
    }
  }
  
  # Simulate future regressors, record means and variances
  means <- variances <- rep(0, hmax)
  xsim <- matrix(0, p+1, hmax)
  xsim[1, ] <- 1
  xsim[2:(p+1), 1] <- rev(y[(length(y)-p+1):length(y)])  
  
  for (hh in 1:hmax){
    means[hh] <- t(xsim[, hh]) %*% betas[, fs[hh]]
    variances[hh] <- vs[fs[hh]]
    if (hh < hmax) xsim[-1, hh+1] <- replace.reg(xsim[-1,hh], rnorm(1, mean = means[hh], sd = sqrt(variances[hh])))    
  }
 
  return(list(m = means, s = sqrt(variances)))
}


##############################################
# Some helper functions for handling time series data
##############################################

tdiff <- function(d1, d2, freq){
  y1 <- substr(d1,1,4)
  y2 <- substr(d2,1,4)
  if (freq == 4) aux <- 6 else aux <- 7
  m1 <- substr(d1,6,aux) 
  m2 <- substr(d2,6,aux)
  return((as.numeric(y2)-as.numeric(y1))*freq+as.numeric(m2)-as.numeric(m1))
}

plust <- function(d1, t, freq){
  if (freq == 4){
   aux1 <- 6
  } else if (freq == 12) {
   aux1 <- 7
  }
  y1 <- as.numeric(substr(d1,1,4))
  p1 <- as.numeric(substr(d1,6,aux1))
  aux <- y1 + (p1-1)/freq + t/freq
  y2 <- floor(aux)
  p2 <- (aux-y2)*freq+1
  if (freq == 4){
	return(paste0(y2,"Q",p2))
  } else if (freq == 12){
	return(paste0(y2,"M", gsub(" ", "0", format(p2, width = 2))))
  }
}

plusq <- function(d1, t) plust(d1, t, 4)
plusm <- function(d1, t) plust(d1, t, 12)
qdiff <- function(d1, d2) tdiff(d1, d2, 4)
mdiff <- function(d1, d2) tdiff(d1, d2, 12)

##############################################
# sel.complete - Select matrix rows with complete entries 
# Input:
# - dat, matrix/data frame
# Output:
# - matrix containing complete rows of dat (if any)
##############################################
sel.complete <- function(dat){
  sel <- apply(dat,1,function(z) !any(is.na(z)))
  dat[sel,] 
}
