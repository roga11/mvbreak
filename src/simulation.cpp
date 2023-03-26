#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// ==============================================================================
//' @title Standard normal errors using box Muller 
//' 
//' @description This function generates uncorrelated standard normal processes using box Muller method.
//' 
//' @param \code{n} Integer determining the length of the process to be simulated
//' @param \code{q}  Integer determining the number of processes to be simulated
//' 
//' @return A (\code{T x q}) matrix of standard normal distributed errors
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
arma::mat randSN(int n, int q){
  double pi = arma::datum::pi;
  arma::mat U1(n, q, arma::fill::randu);
  arma::mat U2(n, q, arma::fill::randu);
  arma::mat eps = trans(trans(sqrt(-2*log(U1))%cos(2*pi*U2)));
  return(eps);
}

// ==============================================================================
//' @title Simulate normally distributed process
//' 
//' @description This function simulates a normally distributed process.
//' 
//' @param mdl_h0 List containing the following DGP parameters
//' \itemize{
//'   \item{\code{n}: }{Length of series.}
//'   \item{\code{mu}: }{A (\code{q x 1}) vector of means.}
//'   \item{\code{sigma}: }{A (\code{q x q}) covariance matrix.}
//'   \item{\code{q}: }{Number of series.}
//'   \item{\code{eps}: }{An optional (\code{T+burnin x q}) matrix with standard normal errors to be used. Errors will be generated if not provided.}
//' }
//' @param burnin Number of simulated observations to remove from beginning. Default is \code{100}.
//' 
//' @return List with simulated series and its DGP parameters.
//' 
//' @export
// [[Rcpp::export]]
List simuNorm(List mdl_h0, int burnin = 0){
  // ----- DGP parameter
  arma::vec mu = mdl_h0["mu"];
  int Tsize = mdl_h0["n"];
  int q = mdl_h0["q"];
  arma::mat eps;
  if(mdl_h0.containsElementNamed("eps") ){
    eps = as<arma::mat>(mdl_h0["eps"]);
  } else {
    eps = randSN(Tsize+burnin, q); 
  }
  // ----- Start simulation
  // pre-define variables
  arma::vec repT(Tsize+burnin, arma::fill::ones);
  // add standard dev & correlations
  arma::mat cov_mat = mdl_h0["sigma"];
  arma::mat C = chol(cov_mat, "lower");
  arma::mat eps_corr = trans(C*trans(eps));
  // simulate process
  arma::mat Y  =  repT*trans(mu) + eps_corr;  
  arma::mat Y_out = Y.submat(burnin,0,Tsize+burnin-1,q-1);
  arma::mat resid_out = eps_corr.submat(burnin,0,Tsize+burnin-1,q-1);
  // ----- Output
  List simuNorm_out = clone(mdl_h0);
  simuNorm_out["y"] = Y_out;
  simuNorm_out["resid"] = resid_out;
  return(simuNorm_out);
}
