#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//' @title Vectorise a matrix row by row
//' 
//' @description Stacks the rows of a matrix into a vector (GAUSS \code{vecr}).
//' 
//' @param x matrix.
//' 
//' @return Numeric vector.
//' 
//' @keywords internal
// [[Rcpp::export]]
arma::vec vecr(arma::mat x){
  arma::mat xt = trans(x);
  arma::vec y = vectorise(xt);
  return(y);
}
  
//' @title Stack a system of equations (C++ core of transf)
//' 
//' @description See \code{transf()}, which calls this.
//' 
//' @param y \code{T} by \code{n} matrix of dependent variables.
//' @param z \code{T} by \code{q} matrix of regressors.
//' @param S selection matrix; ignored when \code{S_NULL} is 1.
//' @param S_NULL 1 to use every regressor in every equation, 0 to apply \code{S}.
//' 
//' @return List with \code{maty} and \code{matz}.
//' 
//' @keywords internal
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
  
//' @title Best single break within a segment
//' 
//' @description Finds the break date in \code{[b1, b2]} that maximises the
//' likelihood of the segment \code{[start, last]} split in two, using the
//' segment likelihoods stored in \code{bigvec}.
//' 
//' @param start first observation of the segment.
//' @param b1,b2 first and last admissible break dates.
//' @param last last observation of the segment.
//' @param bigvec segment likelihoods from \code{mlebigvec()}.
//' @param bigt number of observations.
//' 
//' @return List with \code{loglikmin} (minus the maximised log-likelihood) and
//' \code{bdate} (the break date).
//' 
//' @keywords internal
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

//' @title Best single break within a segment, given the coefficients
//' 
//' @description As \code{parti()}, but the likelihood of each half is computed
//' from the residuals in \code{bigvec2} under given regime coefficients, as
//' in the re-partitioning step of the restricted estimation.
//' 
//' @inheritParams parti
//' @param bigvec2 residuals under the two regimes' coefficients, \code{2T} by \code{n} (from \code{residseg()}).
//' @param n number of equations.
//' 
//' @return List with \code{optmle} (minus the maximised log-likelihood) and \code{bdate}.
//' 
//' @keywords internal
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

//' @title Reshape a matrix row by row
//' 
//' @description Fills an \code{a} by \code{b} matrix row by row with the
//' elements of \code{amat} taken row by row (GAUSS \code{reshape}). Used to
//' turn stacked residuals into a \code{T} by \code{n} matrix.
//' 
//' @param amat matrix whose elements are used.
//' @param a,b dimensions of the result.
//' 
//' @return \code{a} by \code{b} matrix.
//' 
//' @keywords internal
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

//' @title Likelihoods of all segments starting at one date
//' 
//' @description Minus the Gaussian quasi-log-likelihood of the segments
//' \code{[start, last]} for every \code{last} from \code{start + h - 1} to
//' the end of the sample, computed with the recursive updating of Section 4.1
//' of Qu and Perron (2007): iterated feasible GLS on the first \code{h}
//' observations, then one observation at a time with the updating formulas
//' of Tobing and McGilchrist (1992). The values are approximations to the
//' exact maximised likelihoods and depend on \code{h}.
//' 
//' @param start first observation of the segments.
//' @param maty dependent variables stacked by observation.
//' @param matz regressors stacked by observation.
//' @param h minimum length of a segment.
//' @param n number of equations.
//' @param last last observation considered.
//' @param maxit largest number of feasible GLS iterations.
//' 
//' @return Vector of length \code{last}; entries before \code{start + h - 1} are zero.
//' 
//' @references Qu, Z. and P. Perron (2007), "Estimating and Testing Structural Changes in Multivariate Regressions", \emph{Econometrica}, 75, 459-502.
//' @references Tobing, H. and C. A. McGilchrist (1992), "Recursive Residuals for Multivariate Regression Models", \emph{Australian Journal of Statistics}, 34, 217-232.
//' 
//' @keywords internal
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




//' @title Likelihoods of all admissible segments
//' 
//' @description Calls \code{mlef()} for every starting date and stores the
//' results in the triangular layout used by \code{dating_dp()}: the segment
//' starting at \code{i} and ending at \code{j} is at position
//' \code{(i-1)*T + i - (i-1)*i/2 + (j - i)}.
//' 
//' @inheritParams mlef
//' @param bigt number of observations.
//' 
//' @return Vector of length \code{T * (T + 1) / 2}.
//' 
//' @keywords internal
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



//' @title Likelihoods of the segments starting between two dates
//' 
//' @description The part of \code{mlebigvec()} for starting dates
//' \code{start} to \code{end}; used by \code{mlebigvec_par()} to split the
//' work across workers.
//' 
//' @inheritParams mlebigvec
//' @param start,end first and last starting date.
//' 
//' @return Vector holding the entries of \code{mlebigvec()} for those starting dates, in order.
//' 
//' @keywords internal
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


//' @title Dynamic programming over segment likelihoods
//' 
//' @description Finds the partition with the largest likelihood for every number
//' of breaks from 1 to \code{m}, given the vector of segment likelihoods produced
//' by \code{mlebigvec()}. Section 4.1 of Qu and Perron (2007).
//' 
//' @param bigvec vector of minus log-likelihoods of all admissible segments, as
//' returned by \code{mlebigvec()} or \code{mlebigvec_par()}.
//' @param h minimum length of a regime.
//' @param m largest number of breaks considered.
//' @param bigt number of observations.
//' 
//' @return List with \code{glob} (maximised log-likelihood for 1 to \code{m}
//' breaks), \code{datevec} (upper-triangular matrix whose column \code{j} holds
//' the break dates when \code{j} breaks are allowed) and \code{bigvec}.
//' 
//' @keywords internal
// [[Rcpp::export]]
List dating_dp(arma::vec bigvec, int h, int m, int bigt){
  arma::mat datevec(m, m, arma::fill::zeros); 
  arma::mat optdat(bigt, m, arma::fill::zeros);
  arma::mat optmle(bigt, m, arma::fill::zeros);     
  arma::vec dvec(bigt, arma::fill::zeros);
  arma::vec glob(m, arma::fill::zeros);
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

//' @title Global break dates by quasi-maximum likelihood
//' 
//' @description Computes the likelihood of every admissible segment with
//' \code{mlebigvec()} and runs the dynamic programme \code{dating_dp()}.
//' Section 4.1 of Qu and Perron (2007). To compute the segment likelihoods in
//' parallel, call \code{mlebigvec_par()} and \code{dating_dp()} directly.
//' 
//' @param y dependent variables stacked by observation (\code{maty} from \code{transf()}).
//' @param z regressors stacked the same way (\code{matz} from \code{transf()}).
//' @param n number of equations.
//' @param h minimum length of a regime.
//' @param m largest number of breaks considered.
//' @param bigt number of observations.
//' 
//' @return See \code{dating_dp()}.
//' 
//' @references Qu, Z. and P. Perron (2007), "Estimating and Testing Structural Changes in Multivariate Regressions", \emph{Econometrica}, 75, 459-502.
//' 
// [[Rcpp::export]]
List dating_mle(arma::vec y, arma::mat z, int n, int h, int m, int bigt){
  arma::vec bigvec = mlebigvec(y, z, h, n, bigt);
  return(dating_dp(bigvec, h, m, bigt));
}


//' @title Dynamic programming under given regime coefficients
//' 
//' @description The re-partitioning step of the restricted estimation
//' (Section 4.2 of Qu and Perron 2007, step 3): given the residuals of every
//' regime's coefficients over the whole sample, finds the partition that
//' maximises the likelihood for 1 to \code{m} breaks.
//' 
//' @param bigvec2 residuals from \code{residseg()}, \code{T * (m + 1)} by \code{n}.
//' @param h minimum length of a regime.
//' @param m number of breaks.
//' @param n number of equations.
//' @param bigt number of observations.
//' 
//' @return List with \code{glob} (maximised log-likelihood for 1 to \code{m}
//' breaks) and \code{datevec} (break dates, column \code{j} for \code{j} breaks).
//' 
//' @references Qu, Z. and P. Perron (2007), "Estimating and Testing Structural Changes in Multivariate Regressions", \emph{Econometrica}, 75, 459-502.
//' 
//' @keywords internal
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




// Residuals stacked by observation, as a T by n matrix (row t holds the n equations).
static arma::mat resid_matrix(const arma::vec& res, int n){
  return arma::reshape(res, n, res.n_elem / n).t();
}

// The sums X' Omega^{-1} X and X' Omega^{-1} y over the observations, where Omega is
// block diagonal with isig[k] (the inverse covariance of regime k) on the observations
// of that regime; seg holds 0, the break dates and T. The sums are formed observation by
// observation, so no T n by T n matrix is built.
static void gls_sums(const arma::mat& X, const arma::vec& y, int n, const arma::vec& seg,
                     const std::vector<arma::mat>& isig, arma::mat& XtX, arma::vec& Xty){
  XtX.zeros(X.n_cols, X.n_cols); Xty.zeros(X.n_cols);
  for (unsigned int k = 1; k < seg.n_elem; k++){
    for (int t = (int) seg(k - 1); t < (int) seg(k); t++){
      arma::mat Xt = X.rows(t * n, t * n + n - 1);
      arma::mat XtS = Xt.t() * isig[k - 1];
      XtX += XtS * Xt; Xty += XtS * y.subvec(t * n, t * n + n - 1);
    }
  }
}

// The regressors of the model with m breaks: block diagonal expansion of matz by regime
// (the C++ counterpart of pzbar()).
static arma::mat regime_regressors(const arma::mat& matz, int n, const arma::vec& seg){
  int m1 = seg.n_elem - 1, p = matz.n_cols;
  arma::mat out(matz.n_rows, m1 * p, arma::fill::zeros);
  for (int k = 0; k < m1; k++){
    int r0 = (int) seg(k) * n, r1 = (int) seg(k + 1) * n - 1;
    out.submat(r0, k * p, r1, (k + 1) * p - 1) = matz.rows(r0, r1);
  }
  return out;
}

//' @title Unrestricted estimates given the break dates
//' 
//' @description Estimates the coefficients and the error covariance matrix of
//' each regime separately by iterated feasible GLS (equation 16 of Qu and
//' Perron 2007 and the iteration below it).
//' 
//' @param maty dependent variables stacked by observation.
//' @param matz regressors stacked by observation.
//' @param n number of equations.
//' @param m number of breaks.
//' @param br the \code{m} break dates.
//' 
//' @return List with \code{beta} (\code{m + 1} by \code{p} matrix, one row per
//' regime) and \code{vv} (the \code{m + 1} covariance matrices stacked, \code{(m + 1) n} by \code{n}).
//' 
//' @references Qu, Z. and P. Perron (2007), "Estimating and Testing Structural Changes in Multivariate Regressions", \emph{Econometrica}, 75, 459-502.
//' 
//' @keywords internal
// [[Rcpp::export]]
List estim(arma::vec maty, arma::mat matz, int n, int m, arma::vec br){
  arma::mat beta(m+1,matz.n_cols, arma::fill::zeros);
  arma::mat vv((m+1)*n,n,arma::fill::zeros);
  arma::vec seg(m+2,arma::fill::zeros);
  if (m > 0) seg.rows(1,m) = br;
  seg(m+1) = (maty.n_rows/n);
  for (int k = 1; k<=(m+1); k++){
    int i = seg(k-1)+1;
    int j = seg(k);
    double len = j - i + 1;
    arma::mat segx = matz.rows(((i-1)*n+1)-1,(j*n)-1);
    arma::vec segy = maty.rows(((i-1)*n+1)-1,(j*n)-1);
    arma::vec b = pinv(trans(segx)*segx)*trans(segx)*segy;
    arma::mat umat = resid_matrix(segy - segx*b, n);
    arma::mat vvar = trans(umat)*umat/len;
    arma::mat vstar = vvar + 1;
    arma::vec seg1 = {0.0, len};
    arma::mat XtX; arma::vec Xty;
    int itr = 1;
    while (max(abs(vectorise(vvar-vstar)))>1e-6 & itr<200){
      vstar = vvar;
      std::vector<arma::mat> isig(1, pinv(vstar));
      gls_sums(segx, segy, n, seg1, isig, XtX, Xty);
      b = pinv(XtX)*Xty;
      umat = resid_matrix(segy - segx*b, n);
      vvar = trans(umat)*umat/len;
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

//' @title Restricted estimates given the break dates
//' 
//' @description Estimates the regime coefficients under the linear
//' restrictions \code{beta = R delta} and the regime covariance matrices by
//' iterated feasible GLS, starting from the unrestricted estimates of
//' \code{estim()}. Which parameters may change is set by \code{brbeta} and
//' \code{brv}; when only the covariance changes the coefficients are common
//' to all regimes, and when only the coefficients change one covariance matrix
//' is estimated from all regimes.
//' 
//' @inheritParams estim
//' @param bigt number of observations.
//' @param R restriction matrix, \code{beta = R delta}.
//' @param brbeta 1 if the regression coefficients may change, 0 otherwise.
//' @param brv 1 if the covariance matrix may change, 0 otherwise.
//' @param maxit largest number of feasible GLS iterations; a warning is given if reached.
//' 
//' @return List with \code{beta} (the \code{(m + 1) p} stacked coefficients) and
//' \code{vv} (the \code{m + 1} covariance matrices stacked).
//' 
//' @references Qu, Z. and P. Perron (2007), "Estimating and Testing Structural Changes in Multivariate Regressions", \emph{Econometrica}, 75, 459-502.
//' 
//' @keywords internal
// [[Rcpp::export]]
List r_estim(arma::vec maty, arma::mat matz,int bigt, int n, int m, arma::vec br, arma::mat R, int brbeta, int brv, int maxit=1000){
  int bigTs = maty.n_rows;
  arma::vec seg(m+2,arma::fill::zeros);
  if (m > 0) seg.rows(1,m) = br;
  seg(m+1) = (bigTs/n);
  // the unrestricted estimates are the starting values
  List estim_out = estim(maty,matz,n,m,br);
  arma::mat beta = estim_out["beta"];
  arma::mat vv = estim_out["vv"];
  arma::vec nbeta;
  arma::mat nvv;
  arma::mat XtX; arma::vec Xty;
  std::vector<arma::mat> isig(m + 1);
  // inverse covariance matrices, one per regime, from the stacked vv
  auto set_isig = [&](){ for (int k = 0; k <= m; k++) isig[k] = pinv(vv.rows(k*n, (k+1)*n-1)); };
  // covariance matrix of each regime from the residuals of the coefficient vector b on X
  auto update_vv = [&](const arma::mat& X, const arma::vec& b){
    for (int k = 1; k<=(m+1); k++){
      int i = seg(k-1)+1;
      int j = seg(k);
      arma::mat umat = resid_matrix(maty.rows(((i-1)*n+1)-1,(j*n)-1) - X.rows(((i-1)*n+1)-1,(j*n)-1)*b, n);
      vv.rows(((k-1)*n+1)-1,(k*n)-1) = trans(umat)*umat/(double)(j-i+1);
    }
  };
  if (brv==1 & brbeta==1){
    // changes in the coefficients and in the covariance matrix
    arma::mat pmatz = regime_regressors(matz, n, seg)*R;
    set_isig();
    gls_sums(pmatz, maty, n, seg, isig, XtX, Xty);
    arma::vec b = pinv(XtX)*Xty;
    arma::vec bstar = b + 10;
    int itr = 1;
    while (max(abs(bstar-b))>=1e-6 & itr<=maxit){
      bstar = b;
      update_vv(pmatz, bstar);
      set_isig();
      gls_sums(pmatz, maty, n, seg, isig, XtX, Xty);
      b = pinv(XtX)*Xty;
      itr = itr + 1;
    }
    nbeta = R*b;
    nvv = vv;
  }else if (brv==1 & brbeta!=1){
    // changes in the covariance matrix only: common coefficients, regime covariances
    set_isig();
    gls_sums(matz, maty, n, seg, isig, XtX, Xty);
    arma::vec b = pinv(XtX)*Xty;
    arma::vec bstar = b + 10;
    int itr = 1;
    while (max(abs(bstar-b))>=1e-6 & itr<=maxit){
      bstar = b;
      update_vv(matz, bstar);
      set_isig();
      gls_sums(matz, maty, n, seg, isig, XtX, Xty);
      b = pinv(XtX)*Xty;
      itr = itr + 1;
    }
    nbeta = b;
    nvv = vv;
    for (int k = 2; k<=(m+1); k++){
      nbeta = join_cols(nbeta,b);
    }
  }else if (brv!=1 & brbeta==1){
    // changes in the coefficients only: one covariance matrix from all regimes
    arma::mat vvar(n,n,arma::fill::zeros);
    for (int k = 1; k<=(m+1); k++){
      int i = seg(k-1)+1;
      int j = seg(k);
      arma::mat tempx = matz.rows(((i-1)*n+1)-1,(j*n)-1);
      arma::vec tempy = maty.rows(((i-1)*n+1)-1,(j*n)-1);
      arma::mat umat = resid_matrix(tempy - tempx*trans(beta.row(k-1)), n);
      vvar = vvar + trans(umat)*umat;
    }
    double bigtd = bigt;
    vvar = vvar/bigtd;
    arma::mat pmatz = regime_regressors(matz, n, seg)*R;
    arma::vec seg1 = {0.0, (double) bigt};
    std::vector<arma::mat> isig1(1, pinv(vvar));
    gls_sums(pmatz, maty, n, seg1, isig1, XtX, Xty);
    arma::vec b = pinv(XtX)*Xty;
    arma::vec bstar = b + 10;
    int itr = 1;
    while ((max(abs(bstar-b))>=1e-6) & (itr<=maxit)){
      bstar = b;
      arma::mat umat = resid_matrix(maty - pmatz*bstar, n);
      vvar = trans(umat)*umat/bigtd;
      isig1[0] = pinv(vvar);
      gls_sums(pmatz, maty, n, seg1, isig1, XtX, Xty);
      b = pinv(XtX)*Xty;
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

















