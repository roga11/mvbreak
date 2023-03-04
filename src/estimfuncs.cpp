#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//' @title convert mat to vec of transpose
//' 
//' @description This procedure 
//' 
//' @param \code{x}
//' 
//' @return vectorised mat 
//' 
//' @references Qu & Perron (2007)
//' @export
// [[Rcpp::export]]
arma::vec vecr(arma::mat x){
  arma::mat xt = trans(x);
  arma::vec y = vectorise(xt);
  return(y);
}
  
  

//' @title Reshape mat 
//' 
//' @description This procedure 
//' 
//' @param \code{amat} 
//' @param \code{a} 
//' @param \code{b} 
//' 
//' @return reshaped matrix
//' 
//' @references Qu & Perron (2007)
//' @export
// [[Rcpp::export]]
arma::mat reshapeg(arma::mat amat, int a, int b){
  arma::mat bmat(a,b,arma::fill::zeros);
  int r = amat.n_rows;
  int c = amat.n_cols;
  if (r==1 & c==1){
    arma::mat matones(a,b,arma::fill::ones);
    bmat = amat*matones;    
  }else if (r*c >=a*b ){
    arma::mat zvec = trans(amat);
    zvec = vectorise(zvec);
    int kk = 0;
    for (int j = 1; j<=a; j++){
      for (int k = 1; k<=b; k++){
        kk = kk + 1;
        bmat(j-1,k-1) = zvec(kk-1,0);            
      }
    }
  }else if (r*c < a*b){
    int dd=ceil((a*b)/(r*c));
    if (dd >= 2){
      for (int jk = 1; jk<=(dd-1); jk++){
        amat = join_cols(amat,amat);
      }
    }
    arma::mat zvec = trans(amat);
    zvec = vectorise(zvec);
    int kk = 0;
    for (int j = 1; j<=a; j++){
      for (int k = 1; k<=b; k++){
        kk = kk + 1;
        bmat(j-1,k-1) = zvec(kk-1,0);            
      }
    }
  }
  return(bmat);
}

//' @title Computing MLE based on a recursive formula
//' 
//' @description This procedure computes MLE from a data set that starts at date 'start' and ends at date 'last'.
//' 
//' @details Note that that the first h-1 elements have no meaning, usually is initialized as zero.
//' 
//' @param \code{start} starting entry of the sample used.
//' @param \code{last} ending date of the last segment considered.
//' @param \code{maty} vector of dependent variables, with the n*(j-1)+i th element corresponding to the jth observation from ith equation.
//' @param \code{matz} matrix of all regressors included in the regression, arranged in the same way as maty.
//' @param \code{h} minimal length of a segment.
//' @param \code{n} number of equations included in the system.
//' 
//' @return It returns a vector of MLE of length last-start+1 (stored for convenience in a vector of length T).
//' 
//' @references Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
//' @export
// [[Rcpp::export]]
arma::vec mlef(int start, arma::mat maty, arma::mat matz, int h, int n, int last, int maxit = 1000){
  double pi = arma::datum::pi;
  arma::vec vecmle(last, arma::fill::zeros);
  
  // initialize the beta and sigma using the first h observations, 
  // iterated until convergence
  int i       = start;
  int j       = start+h-1;
  arma::mat segx    = matz.rows((i-1)*n,j*n-1);
  arma::mat segy    = maty.rows((i-1)*n,j*n-1);
  arma::vec b       = pinv(trans(segx)*segx)*trans(segx)*segy;
  arma::mat res     = segy - segx*b;
  arma::mat umat    = reshapeg(res,res.n_rows/n,n); 
  arma::mat vvar    = trans(umat)*umat/h;
  arma::mat vstar   = vvar + 1;
  arma::mat diagmat(segy.n_rows/n,segy.n_rows/n,arma::fill::eye);
  arma::mat ibigv;
  int itr = 1;
  while ((max(abs(vecr(vvar-vstar)))>1e-6) and (itr<maxit)){
    vstar = vvar;
    ibigv = kron(diagmat, pinv(vstar));
    b     = pinv(trans(segx)*ibigv*segx)*trans(segx)*ibigv*segy; 
    res   = segy - segx*b;
    umat  = reshapeg(res,res.n_rows/n,n); 
    vvar  = trans(umat)*umat/h;
    itr   = itr + 1;
  }
  arma::mat hdiagmat(h,h,arma::fill::eye);
  ibigv   = kron(hdiagmat,pinv(vvar));
  arma::mat ihstar = pinv(trans(segx)*ibigv*segx);
  arma::vec bstar   = b;
  arma::mat ystar   = segy;
  arma::mat zstar   = segx;
  vstar = vvar;
  arma::mat gstar;
  arma::mat tempz;
  arma::vec tempy;
  arma::mat icstar;
  double idub = i;
  while (j<=last){ // should we start at j+1? (currently counting value at j twice...)
    gstar   = ihstar*trans(matz.rows((j-1)*n,j*n-1));
    icstar  = pinv(vstar + matz.rows((j-1)*n,j*n-1)*ihstar*trans(matz.rows((j-1)*n,j*n-1)));
    b       = bstar + gstar*icstar*(maty.rows((j-1)*n,j*n-1) - matz.rows((j-1)*n,j*n-1)*bstar);
    tempz   = join_cols(zstar, matz.rows((j-1)*n,j*n-1));
    tempy   = join_cols(ystar, maty.rows((j-1)*n,j*n-1));
    res     = tempy - tempz*b;
    umat    = reshapeg(res,res.n_rows/n,n); 
    vvar    = trans(umat)*umat/(j-i+1);
    itr     = 1;
    vstar   = vvar+10;
    while ((max(abs(vecr(vvar-vstar)))>1e-6) & (itr<maxit)){
      vstar   = vvar;
      // note that gstar does not need to be updated
      icstar  = pinv(vstar + matz.rows((j-1)*n,j*n-1)*ihstar*trans(matz.rows((j-1)*n,j*n-1)));
      b       = bstar + gstar*icstar*(maty.rows((j-1)*n,j*n-1) - matz.rows((j-1)*n,j*n-1)*bstar);
      res     = tempy - tempz*b;
      umat    = reshapeg(res,res.n_rows/n,n); 
      vvar    = trans(umat)*umat/(j-i+1);
      itr     = itr + 1;
    }
    double jdub = j; 
    double lentmp = ((jdub-idub+1)/2);
    // note that here we are computing the negative of log-likelihood 
    // because we will be minimizing later.
    vecmle(j-1) = (lentmp*n)*(log(2*pi)+1)+lentmp*log(det(vvar)); 
    bstar     = b;
    zstar     = tempz;
    ystar     = tempy;
    vstar     = vvar;
    ihstar    = ihstar - gstar*icstar*trans(gstar);
    j = j + 1;
  }
  return(vecmle);
}




//' @title MLE function
//' 
//' @references Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
//' @export
// [[Rcpp::export]]
arma::vec mlebigvec(arma::mat maty, arma::mat matz, int h, int n, int bigt){
  arma::vec bigvec(bigt*(bigt+1)/2, arma::fill::zeros);
  arma::vec loglik;
  for (int i =1; i<=(bigt-h+1); i ++){
    loglik = mlef(i, maty, matz, h, n, bigt);
    bigvec.rows(((i-1)*bigt+i-(i-1)*i/2)-1,(i*bigt-(i-1)*i/2)-1) = loglik.rows(i-1,bigt-1);
  }  
  return(bigvec);
}

