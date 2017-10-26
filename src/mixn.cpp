#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;
using namespace R;
using namespace stats;

// [[Rcpp::export]]
double auxcrpsC(double m, double s) {
  if (s <= 0) {
    if (s < 0) {
      return ML_NAN;
    }
    return fabs(m);
  }
	return 2 * s * dnorm_0(m / s, 0) + m * (2 * pnorm_0(m / s, 1, 0) - 1);
}

// [[Rcpp::export]]
double crpsmixnC(NumericVector w, NumericVector m, NumericVector s, double y){
  int N = m.size();
  double crps1 = 0.0;
  double crps2 = 0.0;
  double W = 0.0;
  
  if (w.size() != N || s.size() != N) {
    return ML_NAN;
  }
  
  for (int i = 0; i < N; i++) {
    if (w[i] < 0.0 || s[i] < 0.0) {
      return ML_NAN;
    }
    W += w[i];
    crps1 += w[i] * auxcrpsC(y - m[i], s[i]);
    
    double crps3 = 0.5 * w[i] * auxcrpsC(0.0, M_SQRT2 * s[i]);
    double si2 = s[i] * s[i];
    for (int j = 0; j < i; j++) {
      crps3 += w[j] * auxcrpsC(m[i] - m[j], sqrt(si2 + s[j] * s[j]));
    }
    crps2 += w[i] * crps3;
  }
  return (crps1 - crps2 / W) / W;
}

// [[Rcpp::export]]
NumericVector lsmixnC(NumericVector w, NumericVector m,
                      NumericVector s, NumericVector y) {
  int N = m.size();
  int nrow = y.size();
  double W = 0.0;
  NumericVector ls(nrow);
  
  if (w.size() != N || s.size() != N) {
    return ML_NAN;
  }
  
  for (int i = 0; i < N; i++) {
    if (w[i] < 0.0) {
      return ML_NAN;
    }
    W += w[i];
    
    for (int j = 0; j < nrow; j++){
      ls[j] += w[i] * dnorm(y[j], m[i], s[i], 0);
	  }
  }
  
  return log(W) - log(ls);
}

// [[Rcpp::export]]
NumericVector dssmixnC(NumericVector w, NumericVector m,
                       NumericVector s, NumericVector y) {
  int N = m.size();
  int nrow = y.size();
  double W = 0.0;
  double M = 0.0;
  double V = 0.0;
  NumericVector dss(nrow);
  
  if (w.size() != N || s.size() != N) {
    return ML_NAN;
  }
  
  for (int i = 0; i < N; i++) {
    if (w[i] < 0.0 || s[i] < 0.0) {
      return ML_NAN;
    }
    W += w[i];
    M += w[i] * m[i];
    V += w[i] * (s[i] * s[i] + m[i] * m[i]);
  }
  M = M / W;
  V = V / W - M * M;
  
  for (int i = 0; i < nrow; i++) {
    dss[i] = pow(y[i] - M, 2.0) / V + log(V);
  }
  return dss;
}

// [[Rcpp::export]]
NumericVector dmixnC(NumericVector m, NumericVector s, NumericVector y){
  int N = m.size();
  int nrow = y.size();
  NumericVector out(nrow);
  
  for (int j = 0; j < nrow; j++){
	  for (int i = 0; i < N; i++) {
		out[j] += dnorm(y[j], m[i], s[i], 0);
	  }
  }
  
  return out / N;
}  

// [[Rcpp::export]]
NumericVector pmixnC(NumericVector m, NumericVector s, NumericVector y){
  int N = m.size();
  int nrow = y.size();
  NumericVector out(nrow);
  
  for (int j = 0; j < nrow; j++){
	  for (int i = 0; i < N; i++) {
		out[j] += pnorm(y[j], m[i], s[i], 1, 0);
	  }
  }
  
  return out / N;
}  
