#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Locally ordered breaks (Qu and Perron 2007, Section 5): one break in each of two
// groups of equations, the first group's break at k1 and the second's at k2 >= k1,
// with all coefficients changing and a constant error covariance matrix.

// FGLS for a given pair (k1, k2), iterated from the covariance matrix 'vvar' passed in
// (the OLS fit when 'ols' is true) until it changes by less than 'tol'. Fills beta
// (base coefficients followed by the shifts), vvar and the T by n residual matrix.
static void ordered_fgls(const arma::vec& maty, const arma::mat& matz, int n, int bigt,
                         const arma::uvec& grp1, int k1, int k2, bool ols, double tol, int maxit,
                         arma::vec& beta, arma::mat& vvar, arma::mat& umat){
  int q = matz.n_cols, p2 = 2 * q;
  // regressors of observation t: [Z_t, Z_t[, grp1] 1{t > k1}, Z_t[, !grp1] 1{t > k2}]
  arma::mat X(n * bigt, p2, arma::fill::zeros);
  X.cols(0, q - 1) = matz;
  for (int t = 1; t <= bigt; t++){
    for (int c = 0; c < q; c++){
      bool shift = grp1(c) ? (t > k1) : (t > k2);
      if (shift) X.col(q + c).rows((t - 1) * n, t * n - 1) = matz.col(c).rows((t - 1) * n, t * n - 1);
    }
  }
  arma::mat ivv;
  if (ols){
    beta = arma::solve(X.t() * X, X.t() * maty);
  } else {
    ivv = arma::inv_sympd(vvar);
    arma::mat XtV(p2, p2, arma::fill::zeros); arma::vec Xty(p2, arma::fill::zeros);
    for (int t = 0; t < bigt; t++){
      arma::mat Xt = X.rows(t * n, t * n + n - 1);
      XtV += Xt.t() * ivv * Xt; Xty += Xt.t() * ivv * maty.subvec(t * n, t * n + n - 1);
    }
    beta = arma::solve(XtV, Xty);
  }
  arma::vec res = maty - X * beta;
  umat = arma::reshape(res, n, bigt).t();
  vvar = umat.t() * umat / bigt;
  arma::mat vstar = vvar + 1.0;
  int itr = 1;
  while (arma::abs(vvar - vstar).max() > tol && itr < maxit){
    vstar = vvar;
    ivv = arma::inv_sympd(vvar);
    arma::mat XtV(p2, p2, arma::fill::zeros); arma::vec Xty(p2, arma::fill::zeros);
    for (int t = 0; t < bigt; t++){
      arma::mat Xt = X.rows(t * n, t * n + n - 1);
      XtV += Xt.t() * ivv * Xt; Xty += Xt.t() * ivv * maty.subvec(t * n, t * n + n - 1);
    }
    beta = arma::solve(XtV, Xty);
    res = maty - X * beta;
    umat = arma::reshape(res, n, bigt).t();
    vvar = umat.t() * umat / bigt;
    itr++;
  }
  if (itr >= maxit) Rcpp::warning("ordered_grid: FGLS did not converge for the pair (%d, %d).", k1, k2);
}

// ==============================================================================
//' @title Grid search for locally ordered breaks
//'
//' @description Evaluates the Gaussian likelihood of the two-break model of
//' Section 5 of Qu and Perron (2007) at every admissible pair of dates
//' \code{h < k1 <= k2 <= bigt - h}, where the coefficients of the first group
//' of equations change after \code{k1} and those of the second group after
//' \code{k2}, with a common covariance matrix, by iterated FGLS. The search
//' follows the original code: for each \code{k1} the pairs are visited in
//' increasing \code{k2}, the first from OLS and the others from the covariance
//' matrix of the previous pair.
//'
//' @param maty stacked dependent variables, \code{n bigt} by 1 (from \code{transf()}).
//' @param matz stacked regressors, \code{n bigt} by \code{q}.
//' @param n number of equations.
//' @param bigt number of observations.
//' @param grp1 logical vector of length \code{q}: \code{TRUE} for the
//' coefficients of the first group of equations.
//' @param h minimum segment length.
//' @param max_gap largest \code{k2 - k1} allowed for the test.
//' @param tol convergence tolerance of the FGLS iteration.
//' @param maxit largest number of FGLS iterations per pair.
//'
//' @return A list: \code{loglik}, a \code{bigt} by \code{bigt} matrix holding
//' \code{-bigt/2 log det(Sigma)} at \code{[k1, k2]} and \code{NA} elsewhere;
//' \code{best}, the pair with the largest likelihood (\code{k1}, \code{k2},
//' \code{beta} of length \code{2q}: base coefficients then shifts, \code{vvar},
//' \code{res} the \code{bigt} by \code{n} residuals, \code{loglik}); and
//' \code{best_gap}, the same for the pairs with \code{k2 - k1 <= max_gap}.
//'
//' @references Qu, Z. and P. Perron (2007), "Estimating and Testing Structural Changes in Multivariate Regressions", \emph{Econometrica}, 75, 459-502.
//'
//' @keywords internal
// [[Rcpp::export]]
List ordered_grid(arma::vec maty, arma::mat matz, int n, int bigt, arma::uvec grp1, int h,
                  int max_gap, double tol = 1e-8, int maxit = 1000){
  int q = matz.n_cols;
  if ((int) grp1.n_elem != q) stop("ordered_grid: 'grp1' must have one element per column of 'matz'.");
  arma::mat loglik(bigt, bigt); loglik.fill(arma::datum::nan);
  double best = -arma::datum::inf, best_gap = -arma::datum::inf;
  List out_best, out_gap;
  arma::vec beta; arma::mat vvar, umat;
  for (int k1 = h + 1; k1 <= bigt - h; k1++){
    for (int k2 = k1; k2 <= bigt - h; k2++){
      ordered_fgls(maty, matz, n, bigt, grp1, k1, k2, k2 == k1, tol, maxit, beta, vvar, umat);
      double ll = -bigt / 2.0 * std::log(arma::det(vvar));
      loglik(k1 - 1, k2 - 1) = ll;
      if (ll > best){
        best = ll;
        out_best = List::create(Named("k1") = k1, Named("k2") = k2, Named("beta") = beta,
                                Named("vvar") = vvar, Named("res") = umat, Named("loglik") = ll);
      }
      if (k2 - k1 <= max_gap && ll > best_gap){
        best_gap = ll;
        out_gap = List::create(Named("k1") = k1, Named("k2") = k2, Named("beta") = beta,
                               Named("vvar") = vvar, Named("res") = umat, Named("loglik") = ll);
      }
    }
  }
  return List::create(Named("loglik") = loglik, Named("best") = out_best, Named("best_gap") = out_gap);
}

// ==============================================================================
//' @title FGLS estimates of the model without breaks
//'
//' @description Iterated FGLS of \code{y_t = (I_n \%x\% t(z_t)) S beta + u_t}
//' with one covariance matrix, from OLS until the covariance matrix changes by
//' less than \code{tol}, as the original code does for the null model of the
//' ordered-breaks test.
//'
//' @inheritParams ordered_grid
//'
//' @return List with \code{beta}, \code{vvar}, \code{res} (\code{bigt} by
//' \code{n}) and \code{loglik}, \code{-bigt/2 log det(vvar)}.
//'
//' @keywords internal
// [[Rcpp::export]]
List fgls_null(arma::vec maty, arma::mat matz, int n, int bigt, double tol = 1e-8, int maxit = 1000){
  int p = matz.n_cols;
  arma::vec beta = arma::solve(matz.t() * matz, matz.t() * maty);
  arma::vec res = maty - matz * beta;
  arma::mat umat = arma::reshape(res, n, bigt).t();
  arma::mat vvar = umat.t() * umat / bigt, vstar = vvar + 1.0;
  int itr = 1;
  while (arma::abs(vvar - vstar).max() > tol && itr < maxit){
    vstar = vvar;
    arma::mat ivv = arma::inv_sympd(vvar);
    arma::mat XtV(p, p, arma::fill::zeros); arma::vec Xty(p, arma::fill::zeros);
    for (int t = 0; t < bigt; t++){
      arma::mat Xt = matz.rows(t * n, t * n + n - 1);
      XtV += Xt.t() * ivv * Xt; Xty += Xt.t() * ivv * maty.subvec(t * n, t * n + n - 1);
    }
    beta = arma::solve(XtV, Xty);
    res = maty - matz * beta;
    umat = arma::reshape(res, n, bigt).t();
    vvar = umat.t() * umat / bigt;
    itr++;
  }
  if (itr >= maxit) Rcpp::warning("fgls_null: FGLS did not converge.");
  return List::create(Named("beta") = beta, Named("vvar") = vvar, Named("res") = umat,
                      Named("loglik") = -bigt / 2.0 * std::log(arma::det(vvar)));
}

// ==============================================================================
//' @title Simulated argmax of the limiting process of locally ordered breaks
//'
//' @description Simulates the joint distribution of Theorem 4 of Qu and
//' Perron (2007): the argmax over \code{v1 <= v2} of the two-sided process
//' \code{Z(v1, v2)} built from four independent Wiener processes on
//' \code{[-M, M]}, whose drifts and scales are the ratios of the quadratic forms
//' \code{P1j, P2j, P12j} (before the breaks) and \code{P1jp1, P2jp1, P12jp1}
//' (after). Within each region \code{Z} is the sum of a function of \code{v1}
//' and a function of \code{v2}, so the constrained maximum is found with
//' running maxima in one pass over the grid.
//'
//' @param P1j,P2j,P12j,P1jp1,P2jp1,P12jp1 the quadratic forms of Theorem 4.
//' @param M half-width of the grid.
//' @param step grid spacing.
//' @param B number of replications.
//'
//' @return \code{B} by 2 matrix of \code{(v1, v2)} at the maximum.
//'
//' @references Qu, Z. and P. Perron (2007), "Estimating and Testing Structural Changes in Multivariate Regressions", \emph{Econometrica}, 75, 459-502.
//'
//' @keywords internal
// [[Rcpp::export]]
arma::mat sim_ordered_argmax(double P1j, double P2j, double P12j, double P1jp1, double P2jp1, double P12jp1,
                             double M, double step, int B){
  int N = (int) std::round(M / step);           // grid points on each side of zero
  arma::mat out(B, 2);
  double s2 = std::sqrt(P2j / P1j), s2p = std::sqrt(P2jp1 / P1j), s1p = std::sqrt(P1jp1 / P1j);
  double d_neg2 = (P2j + 2 * P12j) / (2 * P1j), d_pos2 = P2jp1 / (2 * P1j), d_pos1 = (P1jp1 + 2 * P12jp1) / (2 * P1j);
  for (int b = 0; b < B; b++){
    // four Wiener processes on (0, M], at the grid points step, 2 step, ..., M
    arma::mat W = arma::cumsum(arma::randn(N, 4) * std::sqrt(step), 0);
    // v1 < 0: a_neg(i) at v1 = -i step, i = 1..N (B1 on the negative side, W col 0)
    // v1 >= 0: a_pos(i) at v1 = i step, i = 0..N (B4, W col 3)
    // v2 <= 0: b_neg(i) at v2 = -i step, i = 0..N (B2, W col 1)
    // v2 > 0:  b_pos(i) at v2 = i step, i = 1..N (B3, W col 2)
    double zmax = 0.0, v1max = 0.0, v2max = 0.0;      // (0, 0) gives Z = 0
    // region B: v1 < 0, v2 > 0, unconstrained
    double a_neg_max = -arma::datum::inf; int ia = 0;
    for (int i = 1; i <= N; i++){
      double a = W(i - 1, 0) - i * step / 2.0;
      if (a > a_neg_max){ a_neg_max = a; ia = i; }
    }
    double b_pos_max = -arma::datum::inf; int ib = 0;
    for (int i = 1; i <= N; i++){
      double bb = s2p * W(i - 1, 2) - d_pos2 * i * step;
      if (bb > b_pos_max){ b_pos_max = bb; ib = i; }
    }
    if (a_neg_max + b_pos_max > zmax){ zmax = a_neg_max + b_pos_max; v1max = -ia * step; v2max = ib * step; }
    // region A: v1 < 0, v2 <= 0, v1 <= v2: scan v2 from -M upward with the running max of a_neg over v1 <= v2
    double run = -arma::datum::inf; int irun = 0;
    for (int i = N; i >= 0; i--){                    // v2 = -i step
      if (i >= 1){
        double a = W(i - 1, 0) - i * step / 2.0;     // v1 = v2 is admissible when v2 < 0
        if (a > run){ run = a; irun = i; }
      }
      double bb = (i == 0) ? 0.0 : s2 * W(i - 1, 1) - d_neg2 * i * step;
      if (run + bb > zmax){ zmax = run + bb; v1max = -irun * step; v2max = -i * step; }
    }
    // region C: v1 >= 0, v2 > 0, v1 <= v2: scan v2 upward with the running max of a_pos over 0 <= v1 <= v2
    run = 0.0; irun = 0;                              // a_pos(0) = 0
    for (int i = 1; i <= N; i++){                    // v2 = i step
      double a = s1p * W(i - 1, 3) - d_pos1 * i * step;
      if (a > run){ run = a; irun = i; }
      double bb = s2p * W(i - 1, 2) - d_pos2 * i * step;
      if (run + bb > zmax){ zmax = run + bb; v1max = irun * step; v2max = i * step; }
    }
    out(b, 0) = v1max; out(b, 1) = v2max;
  }
  return out;
}
