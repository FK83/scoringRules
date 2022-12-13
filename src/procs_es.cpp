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

// function for variogram score when no weights are specified
// [[Rcpp::export]]
double vsC(arma::colvec y, arma::mat dat, double p){
  
  double out = 0;
  double d = dat.n_rows;
  for (int i = 1; i < (d+1); i++) {
    for (int j = i; j < (d+1); j++) {
      double vdat = mean(pow(abs(dat.row(i-1) - dat.row(j-1)), p));
      double vy = pow(abs(y[i-1]-y[j-1]), p);
      out += 2*pow(vy - vdat, 2.0);
    }
  }
  return (out);
}

// function for variogram score when weights w_vs are specified
// [[Rcpp::export]]
double vsC_w_vs(arma::colvec y, arma::mat dat, arma::mat w_vs, double p){
  
  double out = 0;
  double d = dat.n_rows;
  for (int i = 1; i < (d+1); i++) {
    for (int j = i; j < (d+1); j++) {
      double vdat = mean(pow(abs(dat.row(i-1) - dat.row(j-1)), p));
      double vy = pow(abs(y[i-1]-y[j-1]), p);
      out += 2*w_vs(i-1,j-1)*pow(vy - vdat, 2.0);
    }
  }
  return (out);
}

// variogram score kernel
// [[Rcpp::export]]
double vskernelC(arma::colvec x1, arma::colvec x2, arma::mat w_vs, double p){
  
  double out = 0;
  double d = x1.size();
  for (int i = 1; i < (d+1); i++) {
    for (int j = i; j < (d+1); j++) {
      double vx1 = pow(abs(x1[i-1] - x1[j-1]), p);
      double vx2 = pow(abs(x2[i-1] - x2[j-1]), p);
      out += 2*w_vs(i-1,j-1)*pow(vx1 - vx2, 2.0);
    }
  }
  return(out);
}

// function for variogram score when weights w are specified
// [[Rcpp::export]]
double vsC_w(arma::colvec y, arma::mat dat, arma::mat w_vs, arma::colvec w, double p){
  
  double s1 = 0;
  double m = dat.n_cols;
  for (int i = 1; i < (m+1); i++) {
    s1 += w[i-1]*vskernelC(dat.col(i-1), y, w_vs, p);
  }
  
  double s2 = 0;
  for (int i = 1; i < (m+1); i++) {
    for (int j = 1; j < (m+1); j++) {
      s2 += w[i-1]*w[j-1]*vskernelC(dat.col(i-1), dat.col(j-1), w_vs, p);
    }
  }
  double out = s1 - s2 / 2;
  return (out);
}

// complete function kept for now; however, exported R function "es_sample" uses separate component functions 
// for "XX" and "XY" parts of score (see below)
// [[Rcpp::export]]
double energyscoreC(arma::colvec y, arma::mat dat, NumericVector w){
  
  double s1 = 0;
  double m = dat.n_cols;
  for (int i = 1; i < (m+1); i++) {
    s1 += w[i-1]*euclnormC(dat.col(i-1) - y);
  }
  
  double s2 = 0;
  for (int i = 1; i < (m+1); i++) {
    for (int j = i; j < (m+1); j++) {
      s2 += 2*w[i-1]*w[j-1]*euclnormC(dat.col(i-1) - dat.col(j-1));
    }
  }
  
  double out = s1 - s2 / 2;
  return (out);
  
}

// "XX" part of Energy score
// [[Rcpp::export]]
double esC_xx(arma::mat dat, NumericVector w){
  
  double m = dat.n_cols;
  double out = 0;
  for (int i = 1; i < (m+1); i++) {
    for (int j = i; j < (m+1); j++) {
      out += 2*w[i-1]*w[j-1]*euclnormC(dat.col(i-1) - dat.col(j-1));
    }
  }
  
  return (out);
  
}

// "XY" part of Energy score
// [[Rcpp::export]]
double esC_xy(arma::colvec y, arma::mat dat, NumericVector w){
  
  double out = 0;
  double m = dat.n_cols;
  for (int i = 1; i < (m+1); i++) {
    out += w[i-1]*euclnormC(dat.col(i-1) - y);
  }
  
  return (out);
  
}


// complete function kept for now; however, exported R function "mmds_sample" uses separate component functions 
// for "XX" and "XY" parts of score (see below)
// [[Rcpp::export]]
double mmdscoreC(arma::colvec y, arma::mat dat, NumericVector w){
  
  double s1 = 0;
  double m = dat.n_cols;
  for (int i = 1; i < (m+1); i++) {
    s1 += w[i-1]*exp(-0.5*pow(euclnormC(dat.col(i-1) - y), 2.0));
  }
  
  double s2 = 0;
  for (int i = 1; i < (m+1); i++) {
    s2 += pow(w[i-1], 2.0);
    for (int j = (i+1); j < (m+1); j++) {
      s2 += 2*w[i-1]*w[j-1]*exp(-0.5*pow(euclnormC(dat.col(i-1) - dat.col(j-1)), 2.0));
    }
  }
  
  double out = s2 / 2 - s1;
  return (out);
  
}

// "XX" part of MMD score
// [[Rcpp::export]]
double mmdsC_xx(arma::mat dat, NumericVector w){
  
  double m = dat.n_cols;
  double out = 0;
  for (int i = 1; i < (m+1); i++) {
    out += pow(w[i-1], 2.0);
    for (int j = (i+1); j < (m+1); j++) {
      out += 2*w[i-1]*w[j-1]*exp(-0.5*pow(euclnormC(dat.col(i-1) - dat.col(j-1)), 2.0));
    }
  }

  return (out);
  
}

// "XY" part of MMD score
// [[Rcpp::export]]
double mmdsC_xy(arma::colvec y, arma::mat dat, NumericVector w){
  
  double out = 0;
  double m = dat.n_cols;
  for (int i = 1; i < (m+1); i++) {
    out += w[i-1]*exp(-0.5*pow(euclnormC(dat.col(i-1) - y), 2.0));
  }
  
  return (out);
  
}
