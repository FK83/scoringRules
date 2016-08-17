// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <cmath>
#define _USE_MATH_DEFINES
#include <math.h>
using namespace Rcpp;        

// [[Rcpp::export]]
arma::colvec mvndrawC(arma::colvec mu, arma::mat sig) {

double k = mu.size();
arma::colvec aux = as<arma::colvec>(rnorm(k));
arma::mat csig = arma::chol(sig).t();
arma::colvec out = mu + csig*aux;
return(out);

}
 
 // [[Rcpp::export]]
List carterkohn(arma::mat y, arma::mat Z, arma::mat Ht, arma::mat Qt, double m, double p, double t, arma::colvec B0, arma::mat V0) {
 
arma::colvec bp = B0; // reuses memory and avoids extra copy
arma::mat Vp = V0;
arma::colvec btt = B0;	// initialize now s/t that their scope extends beyond loop
arma::mat Vtt = V0;
arma::mat f = arma::zeros(p, Ht.n_cols);
arma::mat invf = f;
arma::colvec cfe = arma::zeros(y.n_rows);
arma::mat bt = arma::zeros(t,m);
arma::mat Vt = arma::zeros(pow(m,2),t);
arma::colvec loglik = arma::zeros(1);
arma::mat R = arma::zeros(p,p);
arma::mat H = arma::zeros(p,m);
 
for (int i = 1; i < (t+1); i++) {
  R = Ht.rows((i-1)*p,i*p-1);
  H = Z.rows((i-1)*p,i*p-1);

  cfe = y.col(i-1) - H*bp;   // conditional forecast error
  f = H*Vp*H.t() + R;    	  // variance of the conditional forecast error
  invf = f.i();		  // invert only once
  loglik = loglik + log(det(f)) + cfe.t()*invf*cfe;

  btt = bp + Vp*H.t()*invf*cfe;
  Vtt = Vp - Vp*H.t()*invf*H*Vp;

  if (i < t){
    bp = btt;
    Vp = Vtt + Qt;
  }

  bt.row(i-1) = btt.t();
  Vt.col(i-1) = arma::vectorise(Vtt);
}

arma::mat bdraw = arma::zeros(t,m);
bdraw.row(t-1) = mvndrawC(btt,Vtt).t();   

for (int i = 1; i < t; i++) {	// Backward recursions
    arma::colvec bf = bdraw.row(t-i).t();
    btt = bt.row(t-i-1).t();
    Vtt = arma::reshape(Vt.col(t-i-1),m,m);
    f = Vtt + Qt;
    invf = f.i();
    cfe = bf - btt;
    arma::colvec bmean = btt + Vtt*invf*cfe;
    arma::mat bvar = Vtt - Vtt*invf*Vtt;
    bdraw.row(t-i-1) = mvndrawC(bmean, bvar).t();
}

return List::create(Named("loglik") = loglik, Named("bdraws") = bdraw.t());
}

// [[Rcpp::export]]
arma::mat drawsigmaC(arma::colvec yts, arma::colvec qs, arma::colvec ms, arma::colvec u2s, arma::colvec Sigtdraw, arma::mat Zs, arma::mat Wdraw, arma::colvec sigma_prmean, arma::mat sigma_prvar){

double t = yts.n_rows; 

arma::colvec cprw = arma::zeros(7);
arma::colvec statedraw = arma::zeros(t);
for (int i = 1; i < (t+1); i++) {
  arma::colvec prw = arma::zeros(7);
  for (int k = 1; k < 8; k++) {
    prw(k-1) = qs(k-1) * (1/sqrt(2*M_PI*u2s(k-1)))*exp(-0.5*((pow(yts(i-1) - Sigtdraw(i-1) - ms(k-1) + 1.2704,2))/u2s(k-1)));
  }
  cprw = arma::cumsum(prw/arma::sum(prw));
  double trand = as<double>(runif(1));
  double imix = 0;
  if (trand < cprw[0]){
    imix = 1;
  } else if (trand < cprw[1]) {
    imix = 2;
  } else if (trand < cprw[2]) {
    imix = 3;
  } else if (trand < cprw[3]) {
    imix = 4;
  } else if (trand < cprw[4]) {
    imix = 5;
  } else if (trand < cprw[5]) {
    imix = 6;
  } else if (trand < cprw[6]) {
    imix = 7;
  }
  statedraw(i-1) = imix;  
}

arma::mat vart = arma::zeros(t, 1);
arma::mat yts1 = arma::zeros(t, 1);
for (int i = 1; i < (t+1); i++) {
  double imix = statedraw(i-1);
  vart(i-1, 0) = u2s(imix-1);
  yts1(i-1, 0) = yts(i-1) - ms(imix-1) + 1.2704;
}

arma::mat Sigtdraw_new = carterkohn(yts1.t(),Zs,vart,Wdraw,1,1,t,sigma_prmean,sigma_prvar)["bdraws"];

return(Sigtdraw_new.t());
}

// [[Rcpp::export]]
arma::colvec drawbetaC(arma::colvec y, arma::mat Z, arma::colvec s2, arma::colvec betapriorm, arma::mat betaprioriv){

double t = y.n_rows;  
arma::mat w = arma::diagmat( arma::ones(t) / s2 );  
arma::mat betapostv = (betaprioriv + Z.t() * w * Z).i();
arma::colvec betapostm = betapostv * (betaprioriv * betapriorm + Z.t() * w * y);

return(mvndrawC(betapostm, betapostv));

}

// [[Rcpp::export]]
arma::mat makeregs_fcC(arma::mat ydat, double p){
  
  double M = ydat.n_cols;
  arma::mat out = arma::zeros(M,M);
  out.eye();
  arma::mat aux = out;
  
  for (int i = 1; i < (p+1); i++) {
    
    arma::mat tmp = ydat.row(ydat.n_rows-i);
    out = arma::join_rows(out, arma::kron(aux, tmp));
    
  }
  
  return(out);
} 

// [[Rcpp::export]]
List getfcsts(arma::colvec beta, arma::colvec Sigt0, arma::mat Wdraw, arma::mat ydat, double nf, double p){
  
  // Initialize parameters
  arma::colvec Sigtfc = Sigt0;
  arma::mat ystar = ydat;
  arma::colvec fcd = arma::zeros(nf);
  arma::colvec fcm = arma::zeros(nf);
  arma::colvec fcv = arma::zeros(nf);
  arma::colvec auxz = arma::zeros(1);
  
  for (int hh = 1; hh < (nf+1); hh++) { 
    
    // Draw Sigt
    Sigtfc = Sigtfc + mvndrawC(auxz, Wdraw);
    
    // Compute variance
    arma::mat Hfc = exp(Sigtfc);
    
    // save fc variance
    fcv(hh-1) = arma::as_scalar(Hfc); 
    
    // make & save fc mean
    arma::mat Zfc = makeregs_fcC(ystar,p);
    arma::colvec mtemp = Zfc*beta;
    fcm(hh-1) = arma::as_scalar(mtemp);
    
    // draw & save realization
    arma::colvec ytemp = mvndrawC(mtemp, Hfc);
    ystar = arma::join_cols(ystar,ytemp.t());
    fcd(hh-1) = arma::as_scalar(ytemp);
    
  }
  
  return List::create(Named("mean") = fcm, Named("variance") = fcv, Named("draw") = fcd); 
  
}

// [[Rcpp::export]]
arma::mat meye(double n){
  arma::mat aux = arma::zeros(n,n);
  aux.eye();
  return(aux);  
}

// [[Rcpp::export]]
arma::mat matmult(arma::mat x, double nt){
arma::mat out =  meye(x.n_rows);
if (nt == 1){
  out = x;
} else if (nt > 1) {
  arma::mat tempmat = x;
    for (int ii = 1; ii < (nt); ii++) {
      tempmat = tempmat * x;
    }
  out = tempmat;
} 
return(out);
}

// [[Rcpp::export]]
List bvarFcstC(arma::mat b, arma::mat sig, arma::mat y, double nf){
double k = b.n_rows;
double p = (b.n_cols - 1)/k;
arma::colvec bigy = arma::zeros(p*k);
// Just to be sure: Select last p rows from input data matrix
y = y.rows(y.n_rows - p, y.n_rows - 1); 
// Define y vector in companion form
for (int ii = 1; ii < (p+1); ii++) {
  bigy.rows((ii-1)*k,ii*k-1) = y.row(p-ii).t();
}
arma::colvec nu = b.col(0);   // intercept vector 
arma::mat a = b.cols(1,b.n_cols-1); // VAR coefficient matrices
if (p > 1){
  arma::colvec tmp = arma::zeros(k*(p-1));
  nu = arma::join_cols(nu,tmp);
  arma::mat tmp2 = meye(k*(p-1));
  arma::mat tmp3 = arma::zeros(k*(p-1),k);
  a = arma::join_cols(a, arma::join_rows(tmp2, tmp3));  
}
arma::mat om0 = sig;
if (p > 1){
  arma::mat tmp4 = arma::zeros(k, k*(p-1));
  arma::mat tmp5 = arma::zeros(k*(p-1), k*p);
  om0 = arma::join_cols(arma::join_rows(sig, tmp4), tmp5);  
}
arma::mat fcv = om0;
arma::mat aux = meye(k*p);
for (int hh = 1; hh < (nf); hh++) {
  aux = aux + matmult(a, hh);
  fcv = fcv + matmult(a, hh)*om0*(matmult(a, hh).t());
}
arma::mat aux2 = aux*nu + matmult(a, nf)*bigy;
return List::create(Named("mean") = aux2.rows(0, k-1), Named("variance") = fcv(arma::span(0,k-1), arma::span(0, k-1)), Named("bigy") = bigy); 
}  

// [[Rcpp::export]]
double drawMultinomC(NumericVector probs) {
  RNGScope scope;		// ensure RNG gets set/reset
  int k = probs.size();
  double u = as<double>(runif(1)); // draw uniform
  double cumprobs = 0;
  double out = 0;
  for(int ii = 0; ii < k; ii++) {
    cumprobs += probs[ii];
    if (u < cumprobs){
      out = ii + 1;
      break;
    }
  }  
  return(out);
}


// [[Rcpp::export]]
List filterMarkovMixtureC(arma::vec p, arma::mat P, arma::mat lnpdat) {
 
  double T = lnpdat.n_rows;
  double m = lnpdat.n_cols;
  arma::mat filprob = arma::zeros(T, m);
  arma::vec lnl = arma::zeros(T);
  arma::vec auxx1 = arma::zeros(T);
  arma::vec auxx2 = arma::zeros(T);
  arma::vec pit1t1 = p;
  arma::vec pitt1 = p;
  arma::vec pitt = p;
  arma::vec simstate = arma::zeros(T);
  
  // forward recursions, see Eq. (11.9) and the following one in Greenberg (2013)
  for (int i = 1; i < (T+1); i++) {
    pitt1 = P.t() * pit1t1;
    arma::rowvec lnpmax = max(lnpdat.row(i-1)) * arma::ones(1, m);
    pitt = pitt1 % (exp(lnpdat.row(i-1) - lnpmax).t());
    double lkl = sum(pitt);
    pitt = pitt/lkl;
    lnl(i-1) = log(lkl) + lnpmax(0);   
    filprob.row(i-1) = pitt.t();
    pit1t1 = pitt;  
  }
  
  
  // backward recursions, see Eq. (11.8) in Greenberg
  simstate(T-1) = drawMultinomC(wrap(filprob.row(T-1).t()));
  for (int i = (T-1); i > 0; i--) {
    arma::vec p1 = (filprob.row(i-1).t()) % P.col(simstate(i)-1);
    p1 = p1/sum(p1);
    simstate(i-1) = drawMultinomC(wrap(p1));
  }

  return List::create(Named("lnl") = lnl, Named("filprob") = filprob, Named("simstate") = simstate);
}