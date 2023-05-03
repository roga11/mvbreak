
#' @title Vectorize 
#' 
#' @description procedure to jointly estimates the coefficient and variance break dates following section 5.1 of Qu and Perron (2007). The code follows est.m of Qu and Perron (2007)
#' 
#' @references Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
#' 
#' @export
vec <- function(a){
  v <- c(a)  
  return(v)
}


#' @title Diagonal partition of matrix with m breaks
#' 
#' @description procedure to construct the diagonal partition of matz with m break at dates b.
#' 
#' @references Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.
#' @references Bai, Jushan & Pierre Perron (2003), "Computation and Analysis of Multiple Structural Change Models," \emph{Journal of Applied Econometrics}, 18, 1-22.
#' @references Qu & Perron (2007) Econometrica
#' @references PYZ (2020) QE
#' 
#' @export
pzbar <- function(matz,m,bb,bigt){
  nt <- nrow(matz)
  q1 <- ncol(matz)
  n <- nt/bigt
  if (m==0){
    stop("m=0, no break is allowed, procedure terminates")
    matzb <- 0
  }else{
    matzb <- matrix(0,nt,(m+1)*q1)
    matzb[1:(bb[1,1]*n),1:q1] <- matz[1:(bb[1,1]*n),]
    i <- 2
    while (i <=m){
      matzb[(bb[i-1,1]*n+1):(bb[i,1]*n),((i-1)*q1+1):(i*q1)] <- matz[(bb[i-1,1]*n+1):(bb[i,1]*n),]
      i <- i + 1 
    }
    matzb[(bb[m,1]*n+1):nt,(m*q1+1):((m+1)*q1)] <- matz[(bb[m,1]*n+1):nt,]
  }
  return(matzb)
}

#' @title Transform variables
#' 
#' @description Procedure that generates the matrix of regressors (matz) and that of the dependent variable (maty)
#' 
#' @param y dependent variable as a matrix, with the jth row containing observations at time j
#' @param z regressors as a matrix, with the jth row containing
#' 
#' @return List of class \code{} (\code{S3} object) with attributes including:
#' \itemize{
#'   \item{\code{maty}: }{dependent variable as a row vector, with the (i-1)*n+jth row containing ith observation for equation j.}
#'   \item{\code{matz}: }{regressors as a matrix}
#' }
#' 
#' @export  
transf <- function(y,z,S=NULL){
  # WRITE IN C++ SINCE LOOPING OVER T
  if (is.null(S)==TRUE){
    transf_outcpp <- transf_cpp(y,z,matrix(0,0,0),1)  
  }else{
    transf_outcpp <- transf_cpp(y,z,S,0)  
  }
  # In  <- diag(ncol(y))
  # maty <- as.matrix(y[1,])
  # if (is.null(S)==FALSE){
  #   matz <- kronecker(In, t(z[1,]))%*%S 
  # }else{
  #   matz <- kronecker(In, t(z[1,]))
  # }
  # for (i in 2:nrow(y)){
  #   maty <- rbind(maty, as.matrix(y[i,]))
  #   
  #   if (is.null(S)==FALSE){
  #     matz_tmp <- kronecker(In, t(z[i,]))%*%S 
  #   }else{
  #     matz_tmp <- kronecker(In, t(z[i,]))
  #   }
  #   matz <- rbind(matz, matz_tmp)
  # } 
  #out <- list(maty = maty, matz = matz)
  return(transf_outcpp)
}



#' @title optimal break partitions for a given segment
#' 
#' @description procedure to obtain an optimal one break partitions for a segment that 
#' starts at date start and ends at date last. It return the optimal break
#' and the associated minus likelihood value (MML).
#' 
#' @param start beginning of the segment considered
#' @param b1 first possible break date
#' @param b2 last possible break date
#' @param last end of segment considered
#' @param bigvec vector of upper-triangular 
#' 
#' @references Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.
#' @references Bai, Jushan & Pierre Perron (2003), "Computation and Analysis of Multiple Structural Change Models," \emph{Journal of Applied Econometrics}, 18, 1-22.
#' @references Qu Perron (2007) Econometrica
#' @references PYZ (2020) QE
#' 
#' @export
parti_R <- function(start,b1,b2,last,bigvec,bigt){
  # procedure to obtain an optimal one break partitions for a segment that
  # starts at date start and ends at date last. It returns the optimal break
  # and the associated SSR.
  
  dvec <- matrix(0, bigt, 1) 
  
  # for a given ending date, the index is the inserting break point , 
  # the cell contains the minus maximized likelihood function (MML)
  ini <- (start-1)*bigt-(start-2)*(start-1)/2+1 # the starting date of the data
  for (j in b1:b2){
    jj <- j-start
    k <- j*bigt-(j-1)*j/2+last-j
    dvec[j,1] <- bigvec[(ini+jj),1] + bigvec[k,1]
    # the first term: the MML corresponding to data from date start to j
    # The second term: the ML corresponding to data from (j+1) and lasts for last-j
  }
  loglikmin <- t(min(dvec[b1:b2,1]))
  minindcdvec <- which.min(dvec[b1:b2,1])
  dx <- (b1-1)+t(minindcdvec)
  return(list(loglikmin = loglikmin, bdate = dx)) 
}



#' @title optimal break partition for one break
#' 
#' @description obtain an optimal one break partitions for a segment that starts at 
#' date start and ends at date last,under the given coefficients.
#' 
#' @details the procedure is called only when number of breaks is one and \code{bigvec} must be a matrix of size \code{2*bigt} by \code{n} (i.e. residuals from \link{residseg}).
#' 
#' @param start beginning of the segment considered
#' @param b1 first possible break date
#' @param b2 last possible break date
#' @param last end of segment considered
#' @param bigvec vector of upper-triangular 
#' @param bigresvec segmented residual matrix
#' @param n number of equations
#' 
#' @return List of class \code{} (\code{S3} object) with attributes including:
#' \itemize{
#'   \item{\code{optmle}: }{minus the value of the log-likelihood function under the optimal partitons.}
#'   \item{\code{bdate}: }{the break date.}
#' }
#' 
#' @references Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.
#' @references Bai, Jushan & Pierre Perron (2003), "Computation and Analysis of Multiple Structural Change Models," \emph{Journal of Applied Econometrics}, 18, 1-22.
#' @references Qu Perron (2007) Econometrica
#' @references PYZ (2020) QE
#' 
#' @export
partim1_R <- function(start,b1,b2,last,bigvec2,bigt,n){
  dvec <- matrix(0,bigt,1)   
  for (j in b1:b2){
    m1 <- t(as.matrix(bigvec2[start:j,]))%*%as.matrix(bigvec2[start:j,])/(j-start+1)
    m2 <- t(as.matrix(bigvec2[(1*bigt+j+1):(1*bigt+last),]))%*%as.matrix(bigvec2[(1*bigt+j+1):(1*bigt+last),])/(last-j+1)
    dvec[j,1] <- ((last-start+1)*n/2)*(log(2*pi)+1)+((j-start+1)/2)*log(det(m1))+((last-j+1)/2)*log(det(m2)) 
  }
  optmle <- t(min(dvec[b1:b2,1]))
  minindcdvec <- which.min(dvec[b1:b2,1])
  dx <- (b1-1)+t(minindcdvec)
  return(list(optmle = optmle, bdate = dx)) 
}




#' @title Dating global break dates using log-likelihood 
#' 
#' @description Find the optimal partition allowing m breaks (unrestricted model)
#' 
#' @param y dependent variable.
#' @param z regressors
#' @param b estimated coefficients.
#' @param q total number of regression coefficients in each regime.
#' @param n number of equations
#' @param m number of breaks
#' 
#' @return \code{bigresvec}: A (\code{T*(m+1) x n}) matrix with residuals for all the possible segments
#' 
#' @references Qu Perron (2007)
#' 
#' @export
residseg <- function(maty,matz,nbeta,q,n,m){
  bigt <- nrow(maty)/n
  bigresvec <- matrix(0,bigt*(m+1),n)
  for (i in 1:(m+1)){
    bigresvec[((i-1)*bigt+1):(i*bigt),] <- reshapeg((maty - matz%*%as.matrix(nbeta[((i-1)*q+1):(i*q),])),bigt,n)
  }
  return(bigresvec)
}


#' @title MLE function
#' 
#' @references Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
#' @export
mlebigvec_par <- function(maty, matz, h, n, bigt, workers){
  #bigvec <- matrix(0,bigt*(bigt+1)/2, 1)
  N_i_f <- floor(bigt/workers)
  N_worker_i <- matrix(rep(N_i_f, workers), workers, 1)
  if (sum(N_worker_i)<bigt){
    N_worker_i[1:(bigt-(N_i_f*workers))] <- N_worker_i[1:(bigt-(N_i_f*workers))] + 1  
  }
  # ----- Begin parallel simulations
  segs <- c(0,cumsum(N_worker_i))
  segs[length(segs)] <- segs[length(segs)] - h + 1
  bigvec_all <- matrix(0,bigt*(bigt+1)/2,1)
  `%dopar%` <- foreach::`%dopar%`
  wi <- NULL 
  bigvec_all <- foreach::foreach(wi = 1:workers, .inorder = TRUE, .packages = "multiMNbreaks") %dopar% {
    start <- segs[wi] +1
    end <- segs[wi+1] 
    tmp <- mlebigvec_forpar(maty, matz, h, n, bigt, start, end)
    tmp
  }
  bigvec <- as.matrix(c(unlist(bigvec_all),matrix(0,h*(h+1)/2-h,1)))
  return(bigvec)
}



#' @title Dating global break dates using log-likelihood
#' 
#' @description Find the optimal partition allowing m breaks, unrestricted model
#' 
#' @param maty vector dependent variables, with the (i-1)*n+j th element corresponding to the ith obseration from the jth equation included in the regression
#' @param matz matrix of independent variables, arranged in the same manner as maty. The number of columns corresponds to the number of basic parameters included in the regression
#' @param h    minimum number of observations in each segment
#' @param m    number of breaks.
#' @param bigt number of observations in the sample
#' 
#' @return List of class \code{} (\code{S3} object) with attributes including:
#' \itemize{
#'   \item{\code{global}: }{The maximized likelhood function under optimal partitons.}
#'   \item{\code{datevec}: }{the estimated break dates.}
#'   \item{\code{bigvec}: }{the big vector containing the MML for all possible segments.}
#' }
#' 
#' @references Qu Perron (2007)
#' 
#' @export
dating_mle_R <- function(maty,matz,n,h,m,bigt,workers){
  datevec <- matrix(0, m, m)    # up-triangular matrix contains the estimated break dates for break numbers from one to m
  optdat  <- matrix(0, bigt, m) # row index corresponds to the ending dates, column index corresponds to the number of breaks permitted before the ending date the cell contains the optimal final break date
  optmle  <- matrix(0, bigt, m) # same as above, the cell contains the MML corresponding to that break structure
  dvec    <- matrix(0, bigt, 1) # the index is the date after which we inserting the break point. The cell contains the corresponding MML
  glob    <- matrix(0, m, 1)    # Global MML when i breaks are permitted
  
  if (workers > 0){
    bigvec  <- mlebigvec_par(maty, matz, h, n, bigt, workers)
  }else{
    bigvec  <- mlebigvec(maty, matz, h, n, bigt)  
  }
  
  # Section that applies the dynamic programming algorithm to look for the optimal breaks.
  # The first case is when m = 1. Here the dynamic programming algorithm is not
  # needed and one call to the parti() procedure is enough.
  
  if (m == 1){
    parti_out <- parti(1,h,bigt-h,bigt,bigvec,bigt)
    datevec[1,1]  <- parti_out$bdate
    glob[1,1]     <- parti_out$loglikmin
  }else{
    # when m > 1, a dynamic programming algorithm is used.
    # The first step is to obtain the optimal one break partitions for all
    # possible ending dates from 2h to T-mh+1.
    # The optimal dates are stored in a vector optdat.
    # The associated MML are stored in a vector optssr.
    for (j1 in (2*h):bigt){
      parti_out <- parti(1,h,j1-h,j1,bigvec,bigt)           
      optmle[j1,1]  <- parti_out$loglikmin
      optdat[j1,1]  <- parti_out$bdate
    }
    glob[1,1]      <- optmle[bigt,1]
    datevec[1,1]   <- optdat[bigt,1]
    
    # When this is done the algorithm looks for optimal 2,3,... breaks
    # partitions. The index used is ib.
    for (ib in 2:m){
      if (ib == m){
        # if we have reached the highest number of breaks considered, only 
        # one segment is considered, that which ends at the last date of the sample.
        jlast <- bigt
        for (jb in (ib*h):(jlast-h)){
          dvec[jb,1]  <- optmle[jb,ib-1] + bigvec[(jb+1)*bigt-jb*(jb+1)/2,1]
        }
        optmle[jlast,ib]      <- min(dvec[(ib*h):(jlast-h),1])
        minindcdvec           <- which.min(dvec[(ib*h):(jlast-h),1])
        optdat[jlast,ib]      <- (ib*h-1) + t(minindcdvec)
      }else{
        # if we have not reached the highest number of breaks considered, 
        # we need to loop over the last date of the segment, between (ib+1)*h and T.
        for (jlast in ((ib+1)*h):bigt){
          for (jb in (ib*h):(jlast-h)){
            dvec[jb,1]   <- optmle[jb,ib-1] + bigvec[jb*bigt-jb*(jb-1)/2+jlast-jb,1] 
          }
          optmle[jlast,ib] <- min(dvec[(ib*h):(jlast-h),1])
          minindcdvec  <- which.min(dvec[(ib*h):(jlast-h),1])
          optdat[jlast,ib] <- (ib*h-1) + t(minindcdvec)
        }
      }
      # At each time we loop the results with ib breaks are retrieved
      # and printed
      datevec[ib,ib] <- optdat[bigt,ib]
      for (i in 1:(ib-1)){
        xx <- ib-i
        datevec[xx,ib] <- optdat[datevec[xx+1,ib],xx]  
      }
      glob[ib,1] <- optmle[bigt,ib]
    }
  }
  glob <- -glob
  return(list(glob = glob,datevec = datevec,bigvec = bigvec))
}

#' @title Dating global break dates using log-likelihood 
#' 
#' @description Find the optimal partition allowing m breaks
#' 
#' @param bigresvec segmented residual matrix
#' @param h    minimum number of observations in each segment
#' @param m    number of breaks.
#' @param n number of equations.
#' @param bigt number of observations in the sample
#' 
#' @return List of class \code{} (\code{S3} object) with attributes including:
#' \itemize{
#'   \item{\code{global}: }{The maximized likelhood function under optimal partitons.}
#'   \item{\code{datevec}: }{the estimated break dates.}
#' }
#' 
#' @references Qu Perron (2007)
#' 
#' @export
dating_m2_R <- function(bigvec2,h,m,n,bigt){
  datevec <- matrix(0, m, m)    # up-triangular matrix contains the estimated break dates for break numbers from one to m
  optdat  <- matrix(0, bigt, m) # row index corresponds to the ending dates, column index corresponds to the number of breaks permitted before the ending date the cell contains the optimal final break date
  optmle  <- matrix(0, bigt, m) # same as above, the cell contains the MML corresponding to that break structure
  dvec    <- matrix(0, bigt, 1) # the index is the date after which we inserting the break point. The cell contains the corresponding MML
  glob    <- matrix(0, m, 1)    # Global MML when i breaks are permitted
  if (m == 1){
    partim1_out <- partim1(1,h,bigt-h,bigt,bigvec2,bigt,n)
    datevec[1,1] <- partim1_out$bdate
    glob[1,1] <-  partim1_out$optmle
  }else{
    # First loop. Looking for the optimal one break partitions for break 
    # dates between h and T-h. j1 is the last date of the segments.
    for (j1 in (2*h):bigt){
      partim1_out <- partim1(1,h,j1-h,j1,as.matrix(bigvec2[1:(2*bigt),]),bigt,n)        
      optmle[j1,1] <- partim1_out$optmle
      optdat[j1,1] <- partim1_out$bdate
    }
    glob[1,1] <- optmle[bigt,1]
    datevec[1,1] <- optdat[bigt,1]
    # When this is done the algorithm looks for optimal 2,3,... breaks
    # partitions. The index used is ib.
    for (ib in 2:m){
      if (ib == m){
        # if we have reached the highest number of breaks considered, only 
        # one segment is considered, that which ends at the last date of 
        # the sample.
        jlast <- bigt
        for (jb in (ib*h):(jlast-h)){
          pmle <- t(as.matrix(bigvec2[(bigt*m+jb+1):(bigt*(m+1)),]))%*%as.matrix(bigvec2[(bigt*m+jb+1):(bigt*(m+1)),])/(bigt-jb)
          dvec[jb,1] <- optmle[jb,ib-1] + ((bigt-jb+1)*n/2)*(log(2*pi)+1)+((bigt-jb+1)/2)*log(det(pmle))
        }
        optmle[jlast,ib] <- t(min(dvec[(ib*h):(jlast-h),1]))
        minindcdvec <- which.min(dvec[(ib*h):(jlast-h),1])
        optdat[jlast,ib] <- (ib*h-1) + t(minindcdvec)
      }else{
        # if we have not reached the highest number of breaks considered, 
        # we need to loop over the last date of the segment, between (ib+1)*h and T. 
        jlast=(ib+1)*h;
        for (jlast in ((ib+1)*h):bigt){
          for (jb in (ib*h):(jlast-h)){
            pmle<- t(as.matrix(bigvec2[(bigt*ib+jb+1):(bigt*ib+jlast),]))%*%as.matrix(bigvec2[(bigt*ib+jb+1):(bigt*ib+jlast),])/(jlast-jb)
            dvec[jb,1] <- optmle[jb,ib-1] + ((jlast-jb+1)*n/2)*(log(2*pi)+1)+((jlast-jb+1)/2)*log(det(pmle))
          }
          optmle[jlast,ib] <- t(min(dvec[(ib*h):(jlast-h),1]))
          minindcdvec <- which.min(dvec[(ib*h):(jlast-h),1])
          optdat[jlast,ib] <- (ib*h-1) + t(minindcdvec)
        }  
      }
      # At each time we loop the results with ib breaks are retrieved
      # and printed
      datevec[ib,ib] <- optdat[bigt,ib]
      for (i in 1:(ib-1)){
        xx <- ib-i
        datevec[xx,ib] <- optdat[datevec[xx+1,ib],xx]  
      }
      glob[ib,1] <- optmle[bigt,ib]
    } 
  }
  glob <- -glob
  return(list(glob = glob,datevec = datevec))
}






#' @title Model estimation (Unrestricted)
#' 
#' @description procedure that estimates the coefficients under the optimal partitions. 
#' 
#' @details See eq. 16 of Qu & Perron (2007) and the iterative procedure discussed below it.
#' 
#' @param maty dependent variable.
#' @param matz independent variable.
#' @param n number of equations.
#' @param m number of breaks.
#' @param br break dates.
#' 
#' @return List of class \code{} (\code{S3} object) with attributes including:
#' \itemize{
#'   \item{\code{beta}: }{the estimate of coefficients}
#'   \item{\code{vv}: }{}
#' }
#' 
#' @export
estim_R <- function(maty,matz,n,m,br){
  beta <- matrix(0,m+1,ncol(matz))
  vv <- matrix(0,(m+1)*n,n)
  seg <- c(0,br,(nrow(maty)/n))
  for (k in 1:(m+1)){
    i     <- seg[k]+1
    j     <- seg[k+1]
    segx  <- as.matrix(matz[((i-1)*n+1):(j*n),])
    segy  <- as.matrix(maty[((i-1)*n+1):(j*n),])
    b     <- pracma::pinv(t(segx)%*%segx)%*%t(segx)%*%segy
    res   <- segy - segx%*%b
    umat  <- reshapeg(res,nrow(res)/n,n)
    vvar  <- t(umat)%*%umat/(j-i+1)
    vstar <- vvar + 1
    itr   <- 1
    while (max(abs(vecr(vvar-vstar)))>1e-6 & itr<200){
      vstar   <- vvar; 
      ibigv   <- kronecker(diag(nrow(segy)/n),pracma::pinv(vstar))
      b       <- pracma::pinv(t(segx)%*%ibigv%*%segx)%*%t(segx)%*%ibigv%*%segy
      res     <- segy - segx%*%b
      umat    <- reshapeg(res,nrow(res)/n,n)
      vvar    <- t(umat)%*%umat/(j-i+1)
	    itr     <- itr + 1
    }
    beta[k,] <- t(b)
    vv[((k-1)*n+1):(k*n),] <- vvar
  }
  return(list(beta = beta, vv = vv)) 
}


#' @title Model estimation (with restriction)
#' 
#' @description Procedure that estimates the model allowing cross regime restrictions.
#' 
#' @details See eq. 16 of Qu & Perron (2007) and the iterative procedure discussed below it.
#' 
#' @note For the current version of the code, only linear restriction is allowed. 
#' The model includes partial structural change and block-partial structural change models.
#' The estimation of based on iterative feasible GLS
#' 
#' @param maty: dependent variable.
#' @param matz: independent variable.
#' @param n: number of equations.
#' @param m: number of breaks.
#' @param bigt: sample size.
#' @param br: break dates.
#' @param R: the restriction, taking the form beta=R*teta, where teta is the vector.
#' @param of basic parameters of the restricted model.
#' @param brv: brv=1 if the variance-covariance matrix is allowed to change brv=0 otherwise.
#' @param brbeta: brbeta=1 if the coefficients are allowed to change brbeta=0 otherwise.
#' 
#' @return List of class \code{} (\code{S3} object) with attributes including:
#' \itemize{
#'   \item{\code{nbeta}: }{the estimate of coefficients imposing restrictions}
#'   \item{\code{nvv}: }{}
#' }
#' 
#' @export
r_estim_R <- function(maty,matz,bigt,n,m,br,R,brbeta,brv,maxit=1000){
  seg <- c(0,br,(nrow(maty)/n))
  # now estimate an unrestricted model
  estim_out <- estim(maty,matz,n,m,br)
  beta <- estim_out$beta
  vv <- estim_out$vv
  if (brv==1 & brbeta==1){
    # Here we consider changes in coefficients and covariance matrix
    pmatz <- pzbar(matz,m,as.matrix(br),bigt)
    pmatz <- pmatz%*%R
    ibigv <- diag(bigt*n)
    for (k in 1:(m+1)){
      i <- seg[k]+1
      j <- seg[k+1]
      ibigv[((i-1)*n+1):(j*n),((i-1)*n+1):(j*n)] <- kronecker(diag(j-i+1), pracma::pinv(vv[((k-1)*n+1):(k*n),])) 
    }
    b <- pracma::pinv(t(pmatz)%*%ibigv%*%pmatz)%*%t(pmatz)%*%ibigv%*%maty
    bstar <- b + 10
    itr <- 1
    while (max(abs(bstar-b))>=1e-6 & itr<=maxit){
      bstar <- b
      # update the covariance matrix
      for (k in 1:(m+1)){
        i <- seg[k]+1
        j <- seg[k+1]
        tempx <- pmatz[((i-1)*n+1):(j*n),]
        tempy <- maty[((i-1)*n+1):(j*n),]
        res <- tempy - tempx%*%bstar
        umat <- reshapeg(res,nrow(res)/n,n)
        vvar <- t(umat)%*%umat/(j-i+1)
        vv[((k-1)*n+1):(k*n),] <- vvar
      }
      # create the block diagonal matrix of the covariances
      for (k in 1:(m+1)){
        i <- seg[k]+1
        j <- seg[k+1]
        ibigv[((i-1)*n+1):(j*n),((i-1)*n+1):(j*n)] <- kronecker(diag(j-i+1), pracma::pinv(vv[((k-1)*n+1):(k*n),]))
      }
      b <- pracma::pinv(t(pmatz)%*%ibigv%*%pmatz)%*%t(pmatz)%*%ibigv%*%maty
      itr <- itr + 1
    }
    nbeta <- R%*%b
    nvv <- vv
  }else if (brv==1 & brbeta!=1){
    # Here we consider changes in the covariance matrix only, 
    # iterative method is used for the optimization, with the unrestricted 
    # estimate as the prior input. The procedure is iterated until convergence.
    # Non-linear optimization algorithm could be used in place of the iterative procedure.
    ibigv <- diag(bigt*n)
    for (k in 1:(m+1)){
      i <- seg[k]+1
      j <- seg[k+1]
      ibigv[((i-1)*n+1):(j*n),((i-1)*n+1):(j*n)] <- kronecker(diag(j-i+1), pracma::pinv(vv[((k-1)*n+1):(k*n),]))
    }
    b <- pracma::pinv(t(matz)%*%ibigv%*%matz)%*%t(matz)%*%ibigv%*%maty
    bstar <- b + 10
    itr <- 1
    while (max(abs(bstar-b))>=1e-6 & itr<=maxit){
      bstar <- b
      # update the variance-covariance matrix
      for (k in 1:(m+1)){
        i <- seg[k]+1
        j <- seg[k+1]
        tempx <- matz[((i-1)*n+1):(j*n),]
        tempy <- maty[((i-1)*n+1):(j*n),]
        res <- tempy - tempx%*%bstar
        umat <- reshapeg(res,nrow(res)/n,n)
        vvar <- t(umat)%*%umat/(j-i+1)
        vv[((k-1)*n+1):(k*n),] <- vvar
      }
      for (k in 1:(m+1)){
        i <- seg[k]+1
        j <- seg[k+1]
        ibigv[((i-1)*n+1):(j*n),((i-1)*n+1):(j*n)] <- kronecker(diag(j-i+1), pracma::pinv(vv[((k-1)*n+1):(k*n),]))
      }
      b <- pracma::pinv(t(matz)%*%ibigv%*%matz)%*%t(matz)%*%ibigv%*%maty
      itr <- itr + 1
    }
    nbeta <- b
    nvv <- vv
    for (k in 2:(m+1)){
      nbeta <- rbind(nbeta,b)
    }
  }else if (brv!=1 & brbeta==1){
    # Here we consider changes in coefficients only
    vvar <- 0
    for (k in 1:(m+1)){
      i <- seg[k]+1
      j <- seg[k+1]
      tempx <- matz[((i-1)*n+1):(j*n),]
      tempy <- maty[((i-1)*n+1):(j*n),]
      res <- tempy - tempx%*%as.matrix(beta[k,])
      umat <- reshapeg(res,nrow(res)/n,n)
      vvar <- vvar + t(umat)%*%umat
    }
    vvar <- vvar/bigt
    ibigv <- kronecker(diag(bigt),pracma::pinv(vvar))
    pmatz <- pzbar(matz,m,as.matrix(br),bigt)
    pmatz <- pmatz%*%R
    b <- pracma::pinv(t(pmatz)%*%ibigv%*%pmatz)%*%t(pmatz)%*%ibigv%*%maty
    # now iterate till convergence
    bstar <- b+10
    itr <- 1
    while (max(abs(vec(b-bstar)))>=1e-6 & itr<=maxit){
      bstar <- b
      #k=1
      res <- maty - pmatz%*%bstar
      umat <- reshapeg(res,nrow(res)/n,n)
      vvar <- t(umat)%*%umat/bigt
      ibigv <- kronecker(diag(bigt),pracma::pinv(vvar))
      b <- pracma::pinv(t(pmatz)%*%ibigv%*%pmatz)%*%t(pmatz)%*%ibigv%*%maty
      itr <- itr + 1
    }
    nbeta <- R%*%b
    nvv <- vvar
    for (k in 2:(m+1)){
      nvv <- rbind(nvv,vvar)
    }
  }
  return(list(beta = nbeta, vv = nvv))
}

#' @title Model estimation
#' 
#' @description This procedure estimates the breaks and the coefficients using the restrictions, 
#' based on the iterative method described in section 4.2 of Qu & Perron (2007).
#' 
#' @param maty dependent variable
#' @param matz independent variable
#' @param n number of equations.
#' @param m number of breaks.
#' @param bigt sample size
#' @param trm trimming parameter 
#' @param R matrix of restrictions
#' @param brbeta Indicator. 1 if coefficients can change, 0 otherwise.
#' @param brv Indicator. 1 if covariance matrix can change, 0 otherwise.
#' 
#' @return List of class \code{} (\code{S3} object) with attributes including:
#' \itemize{
#'   \item{\code{dx}: }{the break dates}
#'   \item{\code{delta}: }{the coefficients}
#' }
#' 
#' @seealso \code{\link{}}
#' 
#' @export
est <- function(maty,matz,n,m,bigt,trm, R = NULL ,brbeta = 1,brv = 1, h = NULL, workers = 0){
  q <- ncol(matz)
  if (is.null(R)){
    R <- diag(q*(m+1))
  }
  if (is.null(h)){
    h <- round(trm*bigt)  
  }
  tstar <- bigt
  
  datemle_out <- dating_mle(maty,matz,n,h,m,bigt,workers)  
  
  #global <- datemle_out$global # not used 
  datevec <- datemle_out$datevec
  #bigvec <- datemle_out$bigvec # not used
  
  # find the optimal partition of  the sample without restrictions 
  # on the coefficients, serves as the starting value of the 
  # iterative procedure
  br <- datevec[,m] #break dates
  restim_out <- r_estim(maty,matz,bigt,n,m,datevec[,m],R,brbeta,brv)
  rbeta <- restim_out$beta
  rvv <- restim_out$vv  # not used 
  
  # now find optimal partitions under the estimated coefficients
  diff <- 0
  while (diff!=-1){
    bigvec2 <- residseg(maty,matz,rbeta,q,n,m)
    datem2_out <- dating_m2(bigvec2,h,m,n,bigt)
    #global <- datem2_out$global  # not used 
    datevec <- datem2_out$datevec
    #bigvec <- datem2_out$bigvec  # not used 
    if (all(datevec[,m]==br)){
      diff <- -1
      dx <- br 
    }else{
      br <- datevec[,m]
      restim_out <- r_estim(maty,matz,bigt,n,m,br,R,brbeta,brv)
      rbeta <- restim_out$beta
      rvv <- restim_out$vv  # not used 
      diff <-diff + 1               
    }
  }
  return(list(dx = as.matrix(dx), rbeta = rbeta, rvv = rvv))
}
