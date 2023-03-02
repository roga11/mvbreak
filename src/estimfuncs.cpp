#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


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
//' @details Note that that the first h-1 elements have no meaning, usually is initialized as zero and the maximum number of iterations allowed is 1000.
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
arma::mat mlef(int start, arma::mat maty, arma::mat matz, int h, int n, int last){
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
  int itr=1;
  // List vstar_ls;
  // arma::mat vstarinpd;
  // arma::mat ibigv;
  // while ((max(abs(vec(vvar-vstar)))>1e-6) & (itr<1000)){
  //   vstar = vvar;
  //   vstar_ls = invpd(vstar);
  //   vstarinpd = as<arma::mat>(vstar_ls["xinv"]);
  //   ibigv = kron(diagmat, vstarinpd);
  //   segz2invpd_ls = invpd(trans(segz)*ibigv*segz);
  //   segz2invpd = as<arma::mat>(segz2invpd_ls["xinv"]);
  //   b       = segz2invpd*trans(segz)*ibigv*segy;
  //   res    = segy - segz*b;
  //   vvar    = trans(res)*res/h;
  //   itr   = itr + 1;
  // }
  // List vvar_ls = invpd(vvar);
  // arma::mat vvarinvpd = vvar_ls["xinv"];
  // arma::mat hdiagmat(h,h,arma::fill::eye);
  // ibigv   = kron(hdiagmat,vvarinvpd);
  // List ihstar_ls = invpd(trans(segz)*ibigv*segz);
  // arma::mat ihstar = ihstar_ls["xinv"];
  // arma::vec bstar   = b;
  // arma::vec ystar   = segy;
  // arma::mat zstar   = segz;
  // arma::mat gstar;
  // arma::mat tempz;
  // arma::vec tempy;
  // arma::mat icstar;
  // List icstar_ls;
  // vstar   = vvar;
  // double idub = i;
  // while (j<last){ // should we start at j+1 (currently counting value at j twice...)
  //   gstar   = trans(ihstar*z.rows((j-1)+1,j));
  //   icstar_ls = invpd(vstar + z.rows((j-1)+1,j)*ihstar*trans(z.rows((j-1)+1,j)));
  //   icstar = as<arma::mat>(icstar_ls["xinv"]);
  //   b   = bstar + gstar*icstar*(y.rows((j-1)+1,j) - z.rows((j-1)+1,j)*bstar);
  //   tempz   = join_cols(zstar, z.rows((j-1)+1,j));
  //   tempy   = join_cols(ystar, y.rows((j-1)+1,j));
  //   res     = tempy - tempz*b;
  //   vvar    = trans(res)*res/(j-i+1);
  //   itr     = 1;
  //   vstar   = vvar + 10;
  //   while ((max(abs(vec(vvar-vstar)))>1e-6) & (itr<1000)){
  //     vstar = vvar;
  //     icstar_ls = invpd(vstar + z.rows((j-1)+1,j)*ihstar*trans(z.rows((j-1)+1,j)));
  //     icstar = as<arma::mat>(icstar_ls["xinv"]);
  //     b   = bstar + gstar*icstar*(y.rows((j-1)+1,j) - z.rows((j-1)+1,j)*bstar);
  //     res     = tempy - tempz*b;
  //     vvar    = trans(res)*res/(j-i+1);
  //     itr   = itr + 1;
  //   }
  //   double jdub = j; 
  //   double lentmp = (jdub-idub+1)/2;
  //   loglik.row(j) = (lentmp)*(log(2*pi)+1)+lentmp*log(vvar);
  //   bstar     = b;
  //   zstar     = tempz;
  //   ystar     = tempy;
  //   vstar     = vvar;
  //   ihstar    = ihstar - gstar*icstar*trans(gstar);
  //   j = j + 1;
  // }
  return(umat);
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
    loglik = mlef(i,maty,matz,h,n,bigt);
    bigvec.rows(((i-1)*bigt+i-(i-1)*i/2)-1,(i*bigt-(i-1)*i/2)-1) = loglik.rows(i-1,bigt-1);
  }  
  return(bigvec);
}

