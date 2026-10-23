#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// ==============================================================================
//' @title Standard normal draws
//' 
//' @description An \code{n} by \code{q} matrix of independent standard normal
//' draws by the Box-Muller transform, from R's uniform random number generator,
//' so that \code{set.seed()} makes the result reproducible.
//' 
//' @param n number of rows.
//' @param q number of columns.
//' 
//' @return \code{n} by \code{q} matrix.
//' 
//' @keywords internal
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
//' @title Simulate a multivariate normal series
//' 
//' @description Draws \code{n} observations of a \code{q}-variate normal series
//' with mean \code{mu} and covariance matrix \code{sigma}.
//' 
//' @param mdl_h0 list with \code{n} (number of observations), \code{q} (number
//' of series), \code{mu} (mean vector) and \code{sigma} (covariance matrix).
//' An optional element \code{eps}, an \code{n + burnin} by \code{q} matrix of
//' standard normal draws, is used instead of fresh draws.
//' @param burnin number of initial observations to discard.
//' 
//' @return The input list with the simulated series added as \code{y} (\code{n} by
//' \code{q}) and the errors as \code{resid}.
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
//' @title Integer sequence
//' 
//' @description The sequence \code{from, from + 1, ..., to} (GAUSS \code{seqa}).
//' 
//' @param from,to first and last value.
//' @param by unused; kept for the GAUSS signature.
//' 
//' @return Numeric vector.
//' 
//' @keywords internal
//' @export
// [[Rcpp::export]]
arma::vec seqC(int from, int to, double by){
  arma::vec seqtmp(to-from+1,arma::fill::ones);
  seqtmp(0) = from;
  seqtmp = cumsum(seqtmp);
  return(seqtmp);
}



// Largest value of the sup-F functional of Bai and Perron (1998) on one
// discretised Wiener path, for every number of breaks from 1 to M. Wq holds the
// n x q matrix of scaled partial sums; h is the minimum segment length. Breaks
// must satisfy h <= T_1, T_j + h <= T_{j+1}, T_M <= n - h. Entry m-1 of the
// result is the m-break maximum, not yet divided by q*m.
arma::vec supF_path(const arma::mat& Wq, int h, int M){
  int n = Wq.n_rows;
  double nd = n;
  arma::mat optssr(n, M, arma::fill::zeros);
  arma::mat bigmat(n, n, arma::fill::zeros);
  arma::vec glb(M, arma::fill::zeros);
  // bigmat(i-1, j-1): contribution of a break at i followed by a break (or the end) at j
  for (int xi = 1; xi<=(n-h); xi++){
    arma::vec seqtmp = seqC((xi+h), n, 1);
    bigmat.submat(xi-1,(xi+h)-1,xi-1,n-1) = trans(sum(pow((((xi/nd)*Wq.rows((xi+h)-1,n-1))-(seqtmp/nd)*Wq.row(xi-1)),2),1)/(((seqtmp/nd)*(xi/nd))%((seqtmp/nd)-(xi/nd))));
  }
  // one break ending at j1: best first break in [h, j1-h]
  for (int j1 = (2*h); j1<=n; j1++){
    optssr(j1-1,0) = as_scalar(max(bigmat.submat(h-1,j1-1,(j1-h)-1,j1-1)));
  }
  glb(0) = optssr(n-1,0);
  // ib breaks ending at jlast: last break in [ib*h, jlast-h]
  for (int ib = 2; ib<=M; ib++){
    for (int jlast = ((ib+1)*h); jlast<=n; jlast++){
      optssr(jlast-1,ib-1) = as_scalar(max(optssr.submat((ib*h)-1,ib-2,(jlast-h)-1,ib-2) + bigmat.submat((ib*h)-1,jlast-1,(jlast-h)-1,jlast-1)));
    }
    glb(ib-1) = optssr(n-1,ib-1);
  }
  return(glb);
}

//' @title Simulated draws of the sup-F functional for 1 to M breaks
//' 
//' @description Simulates the limiting distribution of the sup-LR test of
//' Qu and Perron (2007) (their Corollary 2; the sup-F limit of Bai and Perron
//' 1998) with \code{q} changing parameters. The Wiener process is approximated
//' by scaled partial sums of \code{n} standard normal draws, and on every path
//' the maximum of the functional over all admissible partitions is computed for
//' 1 to \code{M} breaks. Uses R's random number generator.
//' 
//' @param M largest number of breaks.
//' @param q number of parameters that change at each break.
//' @param h minimum segment length, in steps of the discretised process.
//' @param n number of steps used to approximate the Wiener process.
//' @param B number of simulated paths.
//' 
//' @return A \code{B} by \code{M} matrix; column \code{m} holds the \code{B}
//' draws of the \code{m}-break sup-LR statistic (not divided by \code{q*m}).
//' 
//' @references Bai, J. and P. Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes", \emph{Econometrica}, 66, 47-78.
//' @references Qu, Z. and P. Perron (2007), "Estimating and Testing Structural Changes in Multivariate Regressions", \emph{Econometrica}, 75, 459-502.
//' 
//' @keywords internal
// [[Rcpp::export]]
arma::mat sim_supF_draws(int M, int q, int h, int n, int B){
  arma::mat draws(B, M, arma::fill::zeros);
  for (int xb = 0; xb<B; xb++){
    arma::mat e = randSN(n,q);
    arma::mat Wq = pow(n,-0.5)*cumsum(e, 0);
    draws.row(xb) = trans(supF_path(Wq, h, M));
  }
  return(draws);
}
