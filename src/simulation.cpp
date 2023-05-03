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




//' @title critical values
//' 
//' @export
// [[Rcpp::export]]
arma::vec critval(int k, int q, double trm, int n, int B){
  int h = round(n*trm);
  arma::vec nulldist(B, arma::fill::zeros);
  arma::mat Wq_i;
  arma::mat Wq_ip1;
  double nd = n;
  for (int xb = 0; xb<B; xb++){
    arma::mat e = randSN(n,q);
    arma::mat bigmat(k,n,arma::fill::zeros);
    Wq_ip1 = pow(n,-0.5)*sum(e.rows(0,n-1),0);
    double lam_ip1 = nd/nd;
    for (int xi = (h-1); xi<(n-h); xi++){
      double lam_i = (xi+1)/nd;
      Wq_i = pow(n,-0.5)*sum(e.rows(0,xi),0);
      double tstat_tmp = as_scalar(((lam_i*Wq_ip1) - (lam_ip1*Wq_i))*trans((lam_i*Wq_ip1) - (lam_ip1*Wq_i)))/(lam_i*lam_ip1*(lam_ip1-lam_i));
      bigmat(0,xi) = tstat_tmp;
    }
    arma::mat tstat_vec = bigmat.cols(h-1,n-h-1)/(k*q);
    double tstat = tstat_vec.max();
    nulldist(xb) = tstat;
  }
  return(nulldist);
}



// //' @title parti for critical values
// //' 
// //' @export
// // [[Rcpp::export]]
// List function(int start,int b1,int b2,int last,arma::mat bigmat,int n){
//   // procedure to obtain an optimal one break partitions for a segment that
//   // starts at date start and ends at date last. It returns the optimal break
//   // and the associated SSR.
//   
//   arma::vec dvec(n, arma::fill::zeros);
//   
//   // for a given ending date, the index is the inserting break point , 
//   // the cell contains the minus maximized likelihood function (MML)
//   int ini = (start-1)*n-(start-2)*(start-1)/2+1; // the starting date of the data
//     for (int j = b1; j<=b2; j++){
//       int jj = j-start;
//       int k = j*n-(j-1)*j/2+last-j;
//       dvec(j-1) = bigmat.row((ini+jj)-1) + bigmat.row(k-1);
//       // the first term: the MML corresponding to data from date start to j
//       // The second term: the ML corresponding to data from (j+1) and lasts for last-j
//     }
//     loglikmin <- t(min(dvec[b1:b2,1]))
//       minindcdvec <- which.min(dvec[b1:b2,1])
//       dx <- (b1-1)+t(minindcdvec)
//       return(list(loglikmin = loglikmin, bdate = dx)) 
// }



//' @title critical values
//' 
//' @export
// [[Rcpp::export]]
arma::mat critval2(int k, int q, double trm, int n, int B){
  int h = round(n*trm);
  arma::mat e = randSN(n,q);
  arma::mat bigmat(n*(n+1)/2,q,arma::fill::zeros);
  for (int xj = 1; xj<=(n-h+1); xj++){
    arma::mat tmpvec(n,q,arma::fill::zeros);
    for (int xi = (xj+h-2); xi<n; xi++){
      tmpvec.row(xi) = pow(n,-0.5)*sum(e.rows(xj-1,xi),0);
    }
    bigmat.rows(((xj-1)*n+xj-(xj-1)*xj/2)-1,(xj*n-(xj-1)*xj/2)-1) = tmpvec.rows(xj-1,n-1);
  }
  return(bigmat);
}



//' @title critical values
//' 
//' @export
// [[Rcpp::export]]
arma::vec critval3(int m, int q, double trm, int n, int B){
  double h = round(trm*n);
  double nd = n;
  double qd = q;
  double md = m;
  arma::vec nulldist(B,arma::fill::zeros);
  arma::vec glb(m,arma::fill::zeros);
  for (int xb = 0; xb<B; xb++){
    arma::mat e = randSN(n,q);
    arma::mat Wq = pow(n,-0.5)*cumsum(e, 0);
    //arma::vec glb(m,arma::fill::zeros);
    arma::mat optssr(n,m,arma::fill::zeros);
    for (int j1 = ((m+1)*h); j1<=n; j1++){
      arma::vec tmp(n,arma::fill::zeros);
      for (int xt = (m*h); xt<=(j1-h); xt++){
        double xtd = xt;
        arma::mat ftmp = (((nd/nd)*Wq.row(xt-1))-((xtd/nd)*Wq.row(n-1)));
        tmp(xt-1) = as_scalar(ftmp*trans(ftmp))/((xtd/nd)*(nd/nd)*((nd/nd)-(xtd/nd)));
      }
      optssr(j1-h-1,0) = max(tmp.rows(h-1,j1-h-1));
    }
    glb(0) = optssr(n-h-1,0);
    for (int ib = 2; ib <=m; ib++){
      for (int j1 = (((m-ib+1)*h)+h); j1<=(n-(ib-1)*h); j1++){
        double j1d = j1;
        arma::vec tmp(n,arma::fill::zeros);
        
        for (int xt = ((m-ib+1)*h); xt<=(j1-h); xt++){
          double xtd = xt;
          arma::mat ftmp = (((j1d/nd)*Wq.row(xt-1))-((xtd/nd)*Wq.row(j1-1)));
          tmp(xt-1) = optssr(j1-1,ib-2) + as_scalar(ftmp*trans(ftmp))/((xtd/nd)*(j1d/nd)*((j1d/nd)-(xtd/nd)));
        }
        optssr(j1-h-1,ib-1) = max(tmp.rows((m-ib+1)*h-1,j1-h-1));
      }
      glb(ib-1) = optssr(n-(ib*h)-1,ib-1);
    }
    nulldist(xb) = (1/(qd*md))*glb(m-1);
  }
  return(nulldist);
}


//' @title critical values
//' 
//' @export
// [[Rcpp::export]]
arma::vec critval4(int m, int q, double trm, int n, int B){
  double h = round(trm*n);
  double nd = n;
  double qd = q;
  double md = m;
  arma::vec nulldist(B,arma::fill::zeros);
  arma::vec glb(m,arma::fill::zeros);
  for (int xb = 0; xb<B; xb++){
    arma::mat e = randSN(n,q);
    arma::mat Wq = pow(n,-0.5)*cumsum(e, 0);
    //arma::vec glb(m,arma::fill::zeros);
    arma::mat optssr(n,m,arma::fill::zeros);
    for (int xt = (m*h); xt<=(n-h); xt++){
      double xtd = xt;
      arma::mat ftmp = (((nd/nd)*Wq.row(xt-1))-((xtd/nd)*Wq.row(n-1)));
      optssr(xt-1) = as_scalar(ftmp*trans(ftmp))/((xtd/nd)*(nd/nd)*((nd/nd)-(xtd/nd)));
    }
    glb(0) = max(optssr.col(0));
    for (int ib = 2; ib <=m; ib++){
      for (int j1 = (((m-ib+1)*h)+h); j1<=(n-(ib-1)*h); j1++){
        double j1d = j1;
        arma::vec tmp(n,arma::fill::zeros);
        
        for (int xt = ((m-ib+1)*h); xt<=(j1-h); xt++){
          double xtd = xt;
          arma::mat ftmp = (((j1d/nd)*Wq.row(xt-1))-((xtd/nd)*Wq.row(j1-1)));
          tmp(xt-1) = optssr(j1-1,ib-2) + as_scalar(ftmp*trans(ftmp))/((xtd/nd)*(j1d/nd)*((j1d/nd)-(xtd/nd)));
        }
        optssr(j1-h-1,ib-1) = max(tmp.rows((m-ib+1)*h-1,j1-h-1));
      }
      glb(ib-1) = max(optssr.col(ib-1));
    }
    nulldist(xb) = (1/(qd*md))*glb(m-1);
  }
  return(nulldist);
}





//' @title critical values
//' 
//' @export
// [[Rcpp::export]]
arma::vec seqC(int from, int to, double by){
  arma::vec seqtmp(to-from+1,arma::fill::ones);
  seqtmp(0) = from;
  seqtmp = cumsum(seqtmp);
  return(seqtmp);
}



//' @title critical values
//' 
//' @export
// [[Rcpp::export]]
arma::vec bp_supF_critval(int m, int q, double trm, int n, int B){
  double h = round(trm*n);
  double nd = n;
  double qd = q;
  double md = m;
  arma::vec nulldist(B,arma::fill::zeros);
  for (int xb = 0; xb<B; xb++){
    arma::mat e = randSN(n,q);
    arma::mat Wq = pow(n,-0.5)*cumsum(e, 0);
    arma::vec glb(m,arma::fill::zeros);
    arma::mat optssr(n,m,arma::fill::zeros);
    arma::mat bigmat(n,n,arma::fill::zeros);
    for (int xi = 1; xi<=(n-h); xi++){
      arma::vec seqtmp = seqC((xi+h),n,1);
      bigmat.submat(xi-1,(xi+h)-1,xi-1,n-1) = trans(sum(pow((((xi/nd)*Wq.rows((xi+h)-1,n-1))-(seqtmp/nd)*Wq.row(xi-1)),2),1)/(((seqtmp/nd)*(xi/nd))%((seqtmp/nd)-(xi/nd))));
    }
    if (m==1){
      optssr.submat((m*h)-1,0,n-1,0) = bigmat.submat((m*h)-1,n-1,n-1,n-1);
      glb(0) = as_scalar(max(optssr.submat((m*h)-1,0,(n-h)-1,0)));  
    }else{
      for (int j1 = (2*h); j1<=n; j1++){
        optssr(j1-1,0) = as_scalar(max(bigmat.submat(h-1,j1-1,(j1-h)-1,j1-1)));
      }
      glb(0) = optssr(n-1,0);
      for (int ib = 2; ib<=m; ib++){
        if (ib==m){
          optssr(n-1,ib-1) = as_scalar(max(optssr.submat((h*m)-1,ib-2,n-1,ib-2) + bigmat.submat((h*m)-1,n-1,n-1,n-1)));
        }else{
          for (int jlast = ((ib+1)*h); jlast<=n; jlast++){
            optssr(jlast-1,ib-1) = as_scalar(max(optssr.submat((ib*h)-1,ib-2,jlast-1,ib-2) + bigmat.submat((ib*h)-1,jlast-1,jlast-1,jlast-1)));
          }
        }
        glb(ib-1) = optssr(n-1,ib-1);
      }
    }
    nulldist(xb) = (1/(qd*md))*glb(m-1);
  }
  return(nulldist);
}



//' @title critical values
//' 
//' @export
// [[Rcpp::export]]
arma::vec bp_supSeq_critval(int l, int q, double trm, int n, int B){
  arma::vec nulldist(B,arma::fill::zeros);
  arma::vec Gx = bp_supF_critval(1, q, trm, n, B);
  if (l>0){
    arma::mat indep_1break(B,l+1,arma::fill::zeros);
    indep_1break.col(0) = Gx;
    for (int xl = 1; xl<(l+1); xl++){
      indep_1break.col(xl) = Gx(arma::randperm(B,B));
    }
    nulldist = max(indep_1break,1);
  }else{
    nulldist = Gx;
  }
  return(nulldist);
}




//' @title critical values
//' 
//' @export
// [[Rcpp::export]]
arma::vec bp_UDmax_critval(int M, int q, double trm, int n, int B){
  arma::vec nulldist(B,arma::fill::zeros);
  arma::vec Gx = bp_supF_critval(1, q, trm, n, B);
  if (M>0){
    arma::mat indep_mbreak(B,M,arma::fill::zeros);
    for (int xm = 1; xm<=M; xm++){
      indep_mbreak.col(xm-1) = bp_supF_critval(xm, q, trm, n, B);
    }
    nulldist = max(indep_mbreak,1);
  }
  return(nulldist);
}















