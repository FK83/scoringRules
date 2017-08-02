## helper function for truncated distributions
calcHess_trunc <- function(z, scale, lb, ub,
                           CDF, PDF, PDFp_over_PDF, G, CRPS) {
  a <- ifelse(is.finite(lb), lb, 0)
  b <- ifelse(is.finite(ub), ub, 0)
  PDFquot_a <- ifelse(is.finite(lb), PDFp_over_PDF[, "lb"], 0)
  PDFquot_b <- ifelse(is.finite(ub), PDFp_over_PDF[, "ub"], 0)
  
  Fxa <- CDF[, "lb"] - CDF[, "z"]
  Fxb <- CDF[, "ub"] - CDF[, "z"]
  Gxa <- G[, "lb"] - G[, "z"]
  Gxb <- G[, "ub"] - G[, "z"]
  
  T_y <- CDF[, "lb"] + CDF[, "ub"] - 2 * CDF[, "z"]
  T_a2 <- (z * Fxa - Gxa + CRPS)
  T_a <- -2 * PDF[, "lb"] * T_a2
  T_b2 <- (z * Fxb - Gxb + CRPS)
  T_b <- 2 * PDF[, "ub"] * T_b2
  
  T_yy <- PDF[, "z"]
  T_ya <- -PDF[, "lb"] * Fxb
  T_yb <- PDF[, "ub"] * Fxa
  #T_ay <- T_ya
  T_aa <- PDF[, "lb"]^2 * (z * (Fxb - Fxa) - CRPS + 4 * T_a2 - a) -
    0.5 * PDFquot_a * T_a
  T_ab <- -PDF[, "lb"] * PDF[, "ub"] * (-CRPS + 2 * T_a2 + 2 * T_b2)
  #T_by <- T_yb
  #T_ba <- T_ab
  T_bb <- PDF[, "ub"]^2 * (z * (Fxa - Fxb) - CRPS + 4 * T_b2 + b) -
    0.5 * PDFquot_b * T_b
  
  d2mu <- 2 / scale * (
    T_yy + T_ya + T_yb +
    T_ya + T_aa + T_ab +
    T_yb + T_ab + T_bb
  )
    
  dmu.dsigma <- dsigma.dmu <- 2 / scale * (
    z * (T_yy + T_ya + T_yb) +
    a * (T_ya + T_aa + T_ab) +
    b * (T_yb + T_ab + T_bb)
  )
    
  d2sigma <- 2 / scale * (
    z^2 * T_yy +
    a^2 * T_aa +
    b^2 * T_bb +
    z * a * 2 * T_ya +
    z * b * 2 * T_yb +
    a * b * 2 * T_ab
  )
  
  Hessian <- cbind(d2mu, d2sigma, dmu.dsigma, dsigma.dmu)
  colnames(Hessia) <- c("d2loc", "d2scale", "dloc.dscale", "dscale.dloc")
  rownames(Hessian) <- NULL
  
  return(Hessian)
}
