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
  
//' @title Transform variables
//' 
//' @description Procedure that generates the matrix of regressors (matz) and that of the dependent variable (maty)
//' 
//' @param y dependent variable as a matrix, with the jth row containing observations at time j
//' @param z regressors as a matrix, with the jth row containing
//' 
//' @return List of class \code{} (\code{S3} object) with attributes including:
//' \itemize{
//'  \item{\code{maty}: }{dependent variable as a row vector, with the (i-1)*n+jth row containing ith observation for equation j.}
//'  \item{\code{matz}: }{regressors as a matrix}
//' }
//' 
//' @export  
// [[Rcpp::export]]
List transf_cpp(arma::mat y, arma::mat z, arma::mat S, int S_NULL){
  int Ts = y.n_rows;
  arma::mat In(y.n_cols, y.n_cols, arma::fill::eye);
  arma::vec maty = trans(y.row(0));
  arma::mat matz;
  arma::mat matz_tmp;
  if (S_NULL==0){
    matz = kron(In, z.row(0))*S; 
  }else{
    matz = kron(In, z.row(0));
  }
  for (int i = 1; i<Ts; i++){
    maty = join_cols(maty, trans(y.row(i)));
    if (S_NULL==0){
      matz_tmp = kron(In, z.row(i))*S; 
    }else{
      matz_tmp = kron(In, z.row(i));
    }
    matz = join_cols(matz, matz_tmp);
  } 
  List out;
  out["maty"] = maty;
  out["matz"] = matz;
  return(out);
}
  
//' @title optimal break partitions for a given segment
//' 
//' @description procedure to obtain an optimal one break partitions for a segment that 
//' starts at date start and ends at date last. It return the optimal break
//' and the associated minus likelihood value (MML).
//' 
//' @param start beginning of the segment considered
//' @param b1 first possible break date
//' @param b2 last possible break date
//' @param last end of segment considered
//' @param bigvec vector of upper-triangular 
//' 
//' @references Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.
//' @references Bai, Jushan & Pierre Perron (2003), "Computation and Analysis of Multiple Structural Change Models," \emph{Journal of Applied Econometrics}, 18, 1-22.
//' @references Qu Perron (2007) Econometrica
//' @references PYZ (2020) QE
//' 
//' @export
// [[Rcpp::export]]
List  parti(int start, int b1, int b2, int last, arma::vec bigvec, int bigt){
  // procedure to obtain an optimal one break partitions for a segment that
  // starts at date start and ends at date last. It returns the optimal break
  // and the associated SSR.
    
  arma::vec dvec(bigt, arma::fill::zeros); 
    
  // for a given ending date, the index is the inserting break point , 
  // the cell contains the minus maximized likelihood function (MML)
  int ini = (start-1)*bigt-(start-2)*(start-1)/2+1; // the starting date of the data
  int jj;
  int k;
  for (int j = b1; j<=b2; j++){
    jj = j-start;
    k = j*bigt-(j-1)*j/2+last-j;
    dvec(j-1) = bigvec((ini+jj)-1) + bigvec(k-1);
  // the first term: the MML corresponding to data from date start to j
  // The second term: the ML corresponding to data from (j+1) and lasts for last-j
  }
  double ssrmin = min(dvec.rows(b1-1,b2-1));
  int minindcdvec = dvec.rows(b1-1,b2-1).index_min();
  int dx = (b1-1) + minindcdvec + 1;
  List output;
  output["loglikmin"] = ssrmin;
  output["bdate"] = dx;
  return(output);
}

//' @title optimal break partition for one break
//' 
//' @description obtain an optimal one break partitions for a segment that starts at 
//' date start and ends at date last,under the given coefficients.
//' 
//' @details the procedure is called only when number of breaks is one and \code{bigvec} must be a matrix of size \code{2*bigt} by \code{n} (i.e. residuals from \link{residseg}).
//' 
//' @param start beginning of the segment considered
//' @param b1 first possible break date
//' @param b2 last possible break date
//' @param last end of segment considered
//' @param bigvec vector of upper-triangular 
//' @param bigresvec segmented residual matrix
//' @param n number of equations
//' 
//' @return List of class \code{} (\code{S3} object) with attributes including:
//' \itemize{
//'   \item{\code{optmle}: }{minus the value of the log-likelihood function under the optimal partitons.}
//'   \item{\code{bdate}: }{the break date.}
//' }
//' 
//' @references Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.
//' @references Bai, Jushan & Pierre Perron (2003), "Computation and Analysis of Multiple Structural Change Models," \emph{Journal of Applied Econometrics}, 18, 1-22.
//' @references Qu Perron (2007) Econometrica
//' @references PYZ (2020) QE
//' 
//' @export
// [[Rcpp::export]]
List partim1(int start, int b1, int b2, int last, arma::mat bigvec2, int bigt, int n){
  double pi = arma::datum::pi;
  arma::vec dvec(bigt,arma::fill::zeros);
  double startd = start;
  double lastd = last;
  double nd = n;
  arma::mat m1;
  arma::mat m2;
  for (int j = b1; j<=b2; j++){
    double jd = j;
    m1 = (trans(bigvec2.rows(start-1,j-1))*bigvec2.rows(start-1,j-1))/(jd-startd+1);
    m2 = (trans(bigvec2.rows((1*bigt+j+1)-1,(1*bigt+last)-1))*bigvec2.rows((1*bigt+j+1)-1,(1*bigt+last)-1))/(lastd-jd+1);
    dvec(j-1) = ((lastd-startd+1)*nd/2)*(log(2*pi)+1)+((jd-startd+1)/2)*log(det(m1))+((lastd-jd+1)/2)*log(det(m2));
  }
  double optmle = min(dvec.rows(b1-1,b2-1));
  int minindcdvec = dvec.rows(b1-1,b2-1).index_min();
  int dx = (b1-1) + minindcdvec + 1;
  List output;
  output["optmle"] = optmle;
  output["bdate"] = dx;
  return(output);
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



//' @title MLE function
//' 
//' @references Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
//' @export
// [[Rcpp::export]]
arma::vec mlebigvec_forpar(arma::mat maty, arma::mat matz, int h, int n, int bigt, int start, int end){
  arma::vec bigvec((end-start+1)*(bigt-start+1)-(end-start+1)*(end-start)/2, arma::fill::zeros);
  arma::vec loglik;
  int endtmp = end-start+1;
  for (int i =1; i<=endtmp; i ++){
    loglik = mlef(i+start-1, maty, matz, h, n, bigt);
    //bigvec.rows(((i-1)*bigt+i-(i-1)*i/2)-1,(i*bigt-(i-1)*i/2)-1) = loglik.rows(i+start-2,bigt-1);
    bigvec.rows(((i-1)*(bigt-start+1)+i-(i-1)*i/2)-1,(i*(bigt-start+1)-(i-1)*i/2)-1) = loglik.rows(i+start-2,bigt-1);
  }  
  
  
  //i=2
  //((i-1)*(bigt-(end-start+1)-1)+i-(i-1)*i/2)-1
  //(i*(bigt-(end-start+1)-1)-(i-1)*i/2)-1
  
  //((i*(bigt-(end-start+1)-1)-(i-1)*i/2)-1)-(((i-1)*(bigt-(end-start+1)-1)+i-(i-1)*i/2)-1)+1
  
  
  return(bigvec);
}


//' @title MLE function
//' 
//' @references Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
//' @export
// [[Rcpp::export]]
arma::vec mlebigvec_forpar2(arma::mat maty, arma::mat matz, int h, int n, int bigt, int start, int end){
  //arma::vec bigvec(bigt*(bigt+1)/2, arma::fill::zeros);
  arma::vec bigvec;
  arma::vec loglik;
  for (int i =start; i<=end; i ++){
    loglik = mlef(i, maty, matz, h, n, bigt);
    //bigvec.rows(((i-1)*bigt+i-(i-1)*i/2)-1,(i*bigt-(i-1)*i/2)-1) = loglik.rows(i-1,bigt-1);
    bigvec = join_cols(bigvec,loglik.rows(i-1,bigt-1));
  }  
  return(bigvec);
}

//' @title Compute global break dates using log-likelihood 
//' 
//' @description This is the main procedure which calculates the break points that globally maximize the loglikelihood function. It returns optimal dates and associated log likelihood for all numbers of breaks less than or equal to m. 
//' 
//' @references Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
//' 
//' @export
// [[Rcpp::export]]
List dating_mle(arma::vec y, arma::mat z, int n, int h, int m, int bigt, int workers){
  Rcpp::Environment multimnbreaks("package:multiMNbreaks");
  Rcpp::Function mlebigvec_par = multimnbreaks["mlebigvec_par"];
  arma::mat datevec(m, m, arma::fill::zeros); 
  arma::mat optdat(bigt, m, arma::fill::zeros);
  arma::mat optmle(bigt, m, arma::fill::zeros);     
  arma::vec dvec(bigt, arma::fill::zeros);
  arma::vec glob(m, arma::fill::zeros);
  arma::vec bigvec;// = mlebigvec(y, z, h, n, bigt);
  
  if (workers > 0){
    bigvec  = as<arma::vec>(mlebigvec_par(y, z, h, n, bigt, workers));
  }else{
    bigvec  = mlebigvec(y, z, h, n, bigt);
  }
  
  // Section that applies the dynamic programming algorithm to look for the optimal breaks.
  // The first case is when m = 1. Here the dynamic programming algorithm is not
  // needed and one call to the parti() procedure is enough.
  if (m == 1){
    List parti_out = parti(1,h,bigt-h,bigt,bigvec,bigt);
    datevec(0,0)  = parti_out["bdate"];
    glob(0)     = parti_out["loglikmin"];
  }else{
    // when m > 1, a dynamic programming algorithm is used.
    // The first step is to obtain the optimal one break partitions for all
    // possible ending dates from 2h to T-mh+1.
    // The optimal dates are stored in a vector optdat.
    // The associated MML are stored in a vector optssr.
    for (int j1 = (2*h); j1<=bigt; j1++){
      List parti_out = parti(1,h,j1-h,j1,bigvec,bigt);
      optmle(j1-1,0)  = parti_out["loglikmin"];
      optdat(j1-1,0)  = parti_out["bdate"];
    }
    datevec(0,0)   = optdat(bigt-1,0);
    glob(0)      = optmle(bigt-1,0);
    // When this is done the algorithm looks for optimal 2,3,... breaks
    // partitions. The index used is ib.
    for (int ib = 2; ib <=m; ib++){
      if (ib == m){
        // if we have reached the highest number of breaks considered, only 
        // one segment is considered, that which ends at the last date of the sample.
        int jlast = bigt;
        for (int jb = (ib*h); jb<=(jlast-h); jb++){
          dvec(jb-1)  = optmle(jb-1,ib-2) + bigvec(((jb+1)*bigt-jb*(jb+1)/2)-1);
        }
        optmle(jlast-1,ib-1)    = min(dvec.rows((ib*h)-1,(jlast-h)-1));
        int  minindcdvec        = dvec.rows((ib*h)-1,(jlast-h)-1).index_min();
        optdat(jlast-1,ib-1)    = (ib*h-1) + minindcdvec + 1;
      }else{
        // if we have not reached the highest number of breaks considered, 
        // we need to loop over the last date of the segment, between (ib+1)*h and T.
        for (int jlast = ((ib+1)*h); jlast<=bigt; jlast++){
          for (int jb = (ib*h); jb<=(jlast-h); jb++){
            dvec(jb-1)   = optmle(jb-1,ib-2) + bigvec((jb*bigt-jb*(jb-1)/2+jlast-jb)-1);
          }
          optmle(jlast-1,ib-1)  = min(dvec.rows((ib*h)-1,(jlast-h)-1));
          int minindcdvec       = dvec.rows((ib*h)-1,(jlast-h)-1).index_min();
          optdat(jlast-1,ib-1)  = (ib*h-1) + minindcdvec + 1;
        }
      }
      // At each time we loop the results with ib breaks are retrieved
      // and printed
      datevec(ib-1,ib-1) = optdat(bigt-1,ib-1);
      for (int i = 1; i<=(ib-1); i++){
        int xx = ib-i;
        datevec(xx-1,ib-1) = optdat(datevec(xx,ib-1)-1, xx-1);
      }
      glob(ib-1) = optmle(bigt-1,ib-1);
    }
  }
  glob = -glob;
  // organize output 
  List output; 
  output["glob"] = glob;
  output["datevec"] = datevec;
  output["bigvec"] = bigvec;
  return(output);
}


//' @title Dating global break dates using log-likelihood 
//' 
//' @description This is the main procedure which calculates the break points that globally maximize the loglikelihood function. It returns optimal dates and associated log likelihood for all numbers of breaks less than or equal to m. 
//' 
//' @param bigresvec segmented residual matrix
//' @param h    minimum number of observations in each segment
//' @param m    number of breaks.
//' @param n number of equations.
//' @param bigt number of observations in the sample
//' 
//' @return List of class \code{} (\code{S3} object) with attributes including:
//' \itemize{
//'   \item{\code{global}: }{The maximized likelhood function under optimal partitons.}
//'   \item{\code{datevec}: }{the estimated break dates.}
//' }
//' 
//' @references Qu Perron (2007)
//' 
//' @export
// [[Rcpp::export]]
List dating_m2(arma::mat bigvec2, int h, int m, int n, int bigt){
  double pi = arma::datum::pi;
  arma::mat datevec(m, m, arma::fill::zeros);
  arma::mat optdat(bigt, m, arma::fill::zeros);
  arma::mat optlr(bigt, m, arma::fill::zeros);                    
  arma::vec dvec(bigt, arma::fill::zeros);
  arma::vec glob(m, arma::fill::zeros);
  double bigtd = bigt;
  double nd = n; 
  if (m==1){
    List dating_out = partim1(1,h,bigt-h,bigt,bigvec2,bigt,n);
    datevec(0,0)    = dating_out["bdate"];
    glob(0)         = dating_out["optmle"];
  }else{
    arma::mat pmle;
    // First loop. Looking for the optimal one break partitions for break 
    // dates between h and T-h. j1 is the last date of the segments.
    for (int j1 = (2*h); j1<=bigt; j1++){
      List dating_out =  partim1(1,h,j1-h,j1,bigvec2.rows(0,(2*bigt)-1),bigt,n);
      optlr(j1-1,0)   = dating_out["optmle"];
      optdat(j1-1,0)  = dating_out["bdate"];
    }                                                                      
    datevec(0,0)  = optdat(bigt-1,0);
    glob(0)       = optlr(bigt-1,0);
    // When this is done the algorithm looks for optimal 2,3,... breaks
    // partitions. The index used is ib.
    for (int ib = 2; ib<=m; ib++){
      if (ib==m){
        // if we have reached the highest number of breaks considered, only 
        // one segment is considered, that which ends at the last date of 
        // the sample.
        int jlast = bigt;
        for (int jb = (ib*h); jb<=(jlast-h); jb++){
          double jbd = jb;
          pmle = (trans(bigvec2.rows(bigt*m+jb,(bigt*(m+1))-1))*bigvec2.rows(bigt*m+jb,(bigt*(m+1))-1))/(bigtd-jbd);
          dvec(jb-1) = optlr(jb-1,ib-2) + ((bigtd-jbd+1)*nd/2)*(log(2*pi)+1)+((bigtd-jbd+1)/2)*log(det(pmle));
        }
        optlr(jlast-1,ib-1)  = min(dvec.rows((ib*h)-1,(jlast-h)-1));
        int minindcdev = dvec.rows((ib*h)-1,(jlast-h)-1).index_min();
        optdat(jlast-1,ib-1) = (ib*h-1) + minindcdev + 1;
      
      }else{
        // if we have not reached the highest number of breaks considered, 
        // we need to loop over the last date of the segment, between (ib+1)*h and T. 
        for (int jlast = ((ib+1)*h); jlast<=bigt; jlast++){
          for (int jb = (ib*h); jb<=(jlast-h); jb++){
            double jbd = jb;
            double jlastd = jlast;
            pmle = (trans(bigvec2.rows(bigt*ib+jb,(bigt*ib+jlast)-1))*bigvec2.rows(bigt*ib+jb,(bigt*ib+jlast)-1))/(jlastd-jbd);
            dvec(jb-1) = optlr(jb-1,ib-2) + ((jlastd-jbd+1)*nd/2)*(log(2*pi)+1)+((jlastd-jbd+1)/2)*log(det(pmle));
          }
          optlr(jlast-1,ib-1) = min(dvec.rows((ib*h)-1,(jlast-h)-1));
          int minindcdev = dvec.rows((ib*h)-1,(jlast-h)-1).index_min();
          optdat(jlast-1,ib-1) = (ib*h-1) + minindcdev + 1;
        }
      }
      datevec(ib-1,ib-1) = optdat(bigt-1,ib-1);
      for (int i = 1; i<=(ib-1); i++){
        int xx = ib-i;
        datevec(xx-1,ib-1) = optdat(datevec(xx,ib-1)-1, xx-1);
      }
      glob(ib-1) = optlr(bigt-1,ib-1);
    }
  }
  glob = -glob;
  // organize output 
  List output; 
  output["glob"] = glob;
  output["datevec"] = datevec;
  return(output);
}




//' @title Model estimation (Unrestricted)
//' 
//' @description procedure that estimates the coefficients under the optimal partitions. 
//' 
//' @details See eq. 16 of Qu & Perron (2007) and the iterative procedure discussed below it.
//' 
//' @param maty dependent variable.
//' @param matz independent variable.
//' @param n number of equations.
//' @param m number of breaks.
//' @param br break dates.
//' 
//' @return List of class \code{} (\code{S3} object) with attributes including:
//' \itemize{
//'   \item{\code{beta}: }{the estimate of coefficients}
//'   \item{\code{vv}: }{}
//' }
//' 
//' @export
// [[Rcpp::export]]
List estim(arma::vec maty, arma::mat matz, int n, int m, arma::vec br){
  Rcpp::Environment multimnbreaks("package:multiMNbreaks");
  Rcpp::Function reshapeg = multimnbreaks["reshapeg"];
  arma::mat beta(m+1,matz.n_cols, arma::fill::zeros);
  arma::mat vv((m+1)*n,n,arma::fill::zeros);
  arma::vec seg(m+2,arma::fill::zeros);
  seg.rows(1,m) = br;
  seg(m+1) = (maty.n_rows/n);
  for (int k = 1; k<=(m+1); k++){
    int i = seg(k-1)+1;
    int j = seg(k);
    double id = i;
    double jd = j;
    arma::mat segx = matz.rows(((i-1)*n+1)-1,(j*n)-1);
    arma::vec segy = maty.rows(((i-1)*n+1)-1,(j*n)-1);
    arma::vec b = pinv(trans(segx)*segx)*trans(segx)*segy;
    arma::mat res = segy - segx*b;
    int bigTstmp = res.n_rows;
    arma::mat  umat = as<arma::mat>(reshapeg(res,bigTstmp/n,n));
    arma::mat vvar = trans(umat)*umat/(jd-id+1);
    arma::mat vstar = vvar + 1;
    int itr = 1;
    while (max(abs(vectorise(vvar-vstar)))>1e-6 & itr<200){
      vstar = vvar; 
      arma::mat ibigsegy(segy.n_rows/n,segy.n_rows/n, arma::fill::eye);
      arma::mat ibigv = kron(ibigsegy,pinv(vstar));
      b = pinv(trans(segx)*ibigv*segx)*trans(segx)*ibigv*segy;
      res = segy - segx*b;
      bigTstmp = res.n_rows;
      umat = as<arma::mat>(reshapeg(res,bigTstmp/n,n));
      vvar = trans(umat)*umat/(jd-id+1);
      itr = itr + 1;
    }
    beta.row(k-1) = trans(b);
    vv.rows(((k-1)*n+1)-1,(k*n)-1) = vvar;
  }
  List output;
  output["beta"] = beta;
  output["vv"] = vv;
  return(output);
}






//' @title Model estimation (with restriction)
//' 
//' @description Procedure that estimates the model allowing cross regime restrictions.
//' 
//' @details See eq. 16 of Qu & Perron (2007) and the iterative procedure discussed below it.
//' 
//' @note For the current version of the code, only linear restriction is allowed. 
//' The model includes partial structural change and block-partial structural change models.
//' The estimation of based on iterative feasible GLS
//' 
//' @param maty: dependent variable.
//' @param matz: independent variable.
//' @param n: number of equations.
//' @param m: number of breaks.
//' @param bigt: sample size.
//' @param br: break dates.
//' @param R: the restriction, taking the form beta=R*teta, where teta is the vector.
//' @param of basic parameters of the restricted model.
//' @param brv: brv=1 if the variance-covariance matrix is allowed to change brv=0 otherwise.
//' @param brbeta: brbeta=1 if the coefficients are allowed to change brbeta=0 otherwise.
//' 
//' @return List of class \code{} (\code{S3} object) with attributes including:
//' \itemize{
//'   \item{\code{nbeta}: }{the estimate of coefficients imposing restrictions}
//'   \item{\code{nvv}: }{}
//' }
//' 
//' @export
// [[Rcpp::export]]
List r_estim(arma::vec maty, arma::mat matz,int bigt, int n, int m, arma::vec br, arma::mat R, int brbeta, int brv, int maxit=1000){
  Rcpp::Environment multimnbreaks("package:multiMNbreaks");
  Rcpp::Function estim = multimnbreaks["estim"];
  Rcpp::Function pzbar = multimnbreaks["pzbar"];
  Rcpp::Function reshapeg = multimnbreaks["reshapeg"];
  int bigTs = maty.n_rows;
  arma::vec seg(m+2,arma::fill::zeros);
  seg.rows(1,m) = br;
  seg(m+1) = (bigTs/n);
  // now estimate an unrestricted model
  List estim_out = estim(maty,matz,n,m,br);
  arma::mat beta = estim_out["beta"];
  arma::mat vv = estim_out["vv"];
  arma::vec nbeta;
  arma::mat nvv;
  if (brv==1 & brbeta==1){
    // Here we consider changes in coefficients and covariance matrix
    arma::mat pmatz = as<arma::mat>(pzbar(matz,m,br,bigt));
    pmatz = pmatz*R;
    arma::mat ibigv(bigt*n,bigt*n,arma::fill::eye);
    for (int k = 1; k<=(m+1); k++){
      int i = seg(k-1)+1;
      int j = seg(k);
      arma::mat diagtmp(j-i+1,j-i+1,arma::fill::eye);
      ibigv.submat(((i-1)*n+1)-1,((i-1)*n+1)-1,(j*n)-1,(j*n)-1) = kron(diagtmp, pinv(vv.rows(((k-1)*n+1)-1,(k*n)-1))); 
    }
    arma::vec b = pinv(trans(pmatz)*ibigv*pmatz)*trans(pmatz)*ibigv*maty;
    arma::vec bstar = b + 10;
    int itr = 1;
    while (max(abs(bstar-b))>=1e-6 & itr<=maxit){
      bstar = b;
      // update the covariance matrix
      for (int k = 1; k<=(m+1); k++){
        int i = seg(k-1)+1;
        int j = seg(k);
        double id = i;
        double jd = j;
        arma::mat tempx = pmatz.rows(((i-1)*n+1)-1,(j*n)-1);
        arma::vec tempy = maty.rows(((i-1)*n+1)-1,(j*n)-1);
        arma::mat res = tempy - tempx*bstar;
        int bigTstmp = res.n_rows;
        arma::mat umat = as<arma::mat>(reshapeg(res,bigTstmp/n,n));
        arma::mat vvar = trans(umat)*umat/(jd-id+1);
        vv.rows(((k-1)*n+1)-1,(k*n)-1) = vvar;
      }
      // create the block diagonal matrix of the covariances
      for (int k = 1; k<=(m+1); k++){
        int i = seg(k-1)+1;
        int j = seg(k);
        arma::mat diagtmp(j-i+1,j-i+1,arma::fill::eye);
        ibigv.submat(((i-1)*n+1)-1,((i-1)*n+1)-1,(j*n)-1,(j*n)-1) = kron(diagtmp, pinv(vv.rows(((k-1)*n+1)-1,(k*n)-1))); 
      }
      b = pinv(trans(pmatz)*ibigv*pmatz)*trans(pmatz)*ibigv*maty;
      itr = itr + 1;
    }
    nbeta = R*b;
    nvv = vv;
  }else if (brv==1 & brbeta!=1){
    // Here we consider changes in the covariance matrix only, 
    // iterative method is used for the optimization, with the unrestricted 
    // estimate as the prior input. The procedure is iterated until convergence.
    // Non-linear optimization algorithm could be used in place of the iterative procedure.
    arma::mat ibigv(bigt*n,bigt*n,arma::fill::eye);
    for (int k = 1; k<=(m+1); k++){
      int i = seg(k-1)+1;
      int j = seg(k);
      arma::mat diagtmp(j-i+1,j-i+1,arma::fill::eye);
      ibigv.submat(((i-1)*n+1)-1,((i-1)*n+1)-1,(j*n)-1,(j*n)-1) = kron(diagtmp, pinv(vv.rows(((k-1)*n+1)-1,(k*n)-1))); 
    }
    arma::vec b = pinv(trans(matz)*ibigv*matz)*trans(matz)*ibigv*maty;
    arma::vec bstar = b + 10;
    int itr = 1;
    while (max(abs(bstar-b))>=1e-6 & itr<=maxit){
      bstar = b;
      // update the variance-covariance matrix
      for (int k = 1; k<=(m+1); k++){
        int i = seg(k-1)+1;
        int j = seg(k);
        double id = i;
        double jd = j;
        arma::mat tempx = matz.rows(((i-1)*n+1)-1,(j*n)-1);
        arma::vec tempy = maty.rows(((i-1)*n+1)-1,(j*n)-1);
        arma::mat res = tempy - tempx*bstar;
        int bigTstmp = res.n_rows;
        arma::mat umat = as<arma::mat>(reshapeg(res,bigTstmp/n,n));
        arma::mat vvar = trans(umat)*umat/(jd-id+1);
        vv.rows(((k-1)*n+1)-1,(k*n)-1) = vvar;
      }
      for (int k = 1; k<=(m+1); k++){
        int i = seg(k-1)+1;
        int j = seg(k);
        arma::mat diagtmp(j-i+1,j-i+1,arma::fill::eye);
        ibigv.submat(((i-1)*n+1)-1,((i-1)*n+1)-1,(j*n)-1,(j*n)-1) = kron(diagtmp, pinv(vv.rows(((k-1)*n+1)-1,(k*n)-1))); 
      }
      b = pinv(trans(matz)*ibigv*matz)*trans(matz)*ibigv*maty;
      itr = itr + 1;
    }
    nbeta = b;
    nvv = vv;
    for (int k = 2; k<=(m+1); k++){
      nbeta = join_cols(nbeta,b);
    }
  }else if (brv!=1 & brbeta==1){
    // Here we consider changes in coefficients only
    arma::mat vvar(n,n,arma::fill::zeros);
    for (int k = 1; k<=(m+1); k++){
      int i = seg(k-1)+1;
      int j = seg(k);
      arma::mat tempx = matz.rows(((i-1)*n+1)-1,(j*n)-1);
      arma::vec tempy = maty.rows(((i-1)*n+1)-1,(j*n)-1);
      arma::mat res = tempy - tempx*trans(beta.row(k-1));
      int bigTstmp = res.n_rows;
      arma::mat umat = as<arma::mat>(reshapeg(res,bigTstmp/n,n));
      vvar = vvar + trans(umat)*umat;
    }
    double bigtd = bigt;
    vvar = vvar/bigtd;
    arma::mat ibigvtmp(bigt,bigt,arma::fill::eye);
    arma::mat ibigv = kron(ibigvtmp,pinv(vvar));
    
    arma::mat pmatz = as<arma::mat>(pzbar(matz,m,br,bigt));
    pmatz = pmatz*R;
    
    arma::vec b = pinv(trans(pmatz)*ibigv*pmatz)*trans(pmatz)*ibigv*maty;
    arma::vec bstar = b + 10;
    int itr = 1;
    while ((max(abs(bstar-b))>=1e-6) & (itr<=maxit)){
      bstar = b;
      arma::mat res = maty - pmatz*bstar;
      int bigTstmp = res.n_rows;
      arma::mat umat = as<arma::mat>(reshapeg(res,bigTstmp/n,n));
      vvar = trans(umat)*umat/bigtd;
      arma::mat ibigvtmp(bigt,bigt,arma::fill::eye);
      arma::mat ibigv = kron(ibigvtmp,pinv(vvar));
      b = pinv(trans(pmatz)*ibigv*pmatz)*trans(pmatz)*ibigv*maty;
      itr = itr + 1;
    }
    nbeta = R*b;
    nvv = vvar;
    for (int k = 2; k<=(m+1); k++){
      nvv = join_cols(nvv,vvar);
    }
  }
  List output;
  output["beta"] = nbeta;
  output["vv"] = nvv;
  return(output);
}

















