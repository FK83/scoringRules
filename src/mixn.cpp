#include <Rcpp.h>
#include <cmath>
#define _USE_MATH_DEFINES
#include <math.h>

using namespace Rcpp;

// [[Rcpp::export]]
double dnormC(double x){
	return (1/sqrt(2*M_PI))*exp(-pow(x,2.0)*0.5);
}

// [[Rcpp::export]]
double pnormC(double x){
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);

    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return 0.5*(1.0 + sign*y);
}                 

// [[Rcpp::export]]
double auxcrpsC(double m, double s){
	return 2*s*dnormC(m/s) + m*(2*pnormC(m/s)-1);
}

// [[Rcpp::export]]
double crpsmixnC(NumericVector w, NumericVector m, NumericVector s, double y){
  int N = m.size();
  double crps = 0;

  for(int i = 0; i < N; i++) {
	crps += (auxcrpsC(y-m[i],s[i])*w[i]);
    for (int j = 0; j < (i+1); j++){
      double tmp = 0.5*auxcrpsC(m[i]-m[j],sqrt(pow(s[i],2.0)+pow(s[j],2.0)))*(w[i]*w[j]);
      crps -= tmp;
      if (i != j){
        crps -= tmp;
      }
    }
  }
  return crps;
}

// [[Rcpp::export]]
NumericVector lsmixnC(NumericVector w, NumericVector m, NumericVector s, NumericVector y){
  int N = m.size();
  int nrow = y.size();
  NumericVector ls(nrow);
  
  for (int j = 0; j < nrow; j++){
	  for(int i = 0; i < N; i++) {
		ls[j] -= (w[i]/s[i])*dnormC((y[j]-m[i])/s[i]);
	  }
  }
  
  return log(ls);
}
