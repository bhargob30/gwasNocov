#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat V_lambda_tau(arma::mat K, arma::vec log_lambda_tau) {
  int R=K.n_rows ;
  return   pow(exp(log_lambda_tau[1]),-1)*(exp(log_lambda_tau[0])*K  + arma:: eye(R,R) )   ;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec beta_tilda(arma::mat K,arma::vec x,arma::vec y, arma::vec log_lambda_tau ){
  arma::mat A= V_lambda_tau(K, log_lambda_tau);
  arma::mat B= arma::inv_sympd(A);
  return( arma::inv_sympd(x.t() * B *x)  * (x.t() * B *y)           );
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double neglr_lambda_tau(arma::vec log_lambda_tau, arma::mat K,arma::vec x,arma::vec y ){
  arma::mat A= V_lambda_tau(K, log_lambda_tau);
  arma::mat B= arma::inv_sympd(A);
  arma::vec b= beta_tilda(K,x,y,log_lambda_tau);
  arma::mat C= 0.5*(arma::log_det_sympd(A)+ log((x.t() * B *x)) + (y- x*b ).t() * B  * (y- x*b ) )  ;
  return(C(0,0));
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec beta_hat_sigma_hat(arma:: mat K, arma::vec x, arma::vec y, arma::vec log_lambda_tau_hat ){
  arma::mat A= V_lambda_tau(K, log_lambda_tau_hat);
  arma::mat B= arma::inv_sympd(A);
  arma::mat C=x.t()*B* x ;
  arma::mat D=arma::inv_sympd(C);
  arma::vec E= D * x.t() * B * y;
  arma::vec out={E(0), pow( D(0, 0),0.5)};
  return out;
  }