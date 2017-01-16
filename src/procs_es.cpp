// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <cmath>
#define _USE_MATH_DEFINES
#include <math.h>
using namespace Rcpp;     

// [[Rcpp::export]]
double euclnormC(arma::colvec x){
  double out = sqrt(sum(square(x)));
  return(out);
}

// [[Rcpp::export]]
double energyscoreC(arma::colvec y, arma::mat dat){
  
  double s1 = 0;
  double m = dat.n_cols;
  for (int i = 1; i < (m+1); i++) {
    s1 += euclnormC(dat.col(i-1) - y);
  }
  
  double s2 = 0;
  for (int i = 1; i < (m+1); i++) {
    for (int j = i; j < (m+1); j++) {
      s2 += 2*euclnormC(dat.col(i-1) - dat.col(j-1));
    }
  }
  
  double out = (s1 / m) - s2 / (2 * pow(m, 2));
  return (out);
  
}