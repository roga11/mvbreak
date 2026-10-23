#' @title Symmetric square root of a positive definite matrix
#'
#' @description Computes \code{M^(1/2)} from the eigendecomposition of \code{M}.
#'
#' @param M symmetric positive definite matrix.
#'
#' @return Symmetric matrix \code{A} with \code{A \%*\% A = M}.
#'
#' @keywords internal
sqrm <- function(M){
  eg <- eigen(M)
  dd <- diag(eg$values)
  ve <- eg$vectors
  sq <- ve%*%sqrt(dd)%*%t(ve)
  return(sq)
}

#' @title Least-squares coefficients
#' 
#' @description Ordinary least squares by QR decomposition, as GAUSS's
#' \code{olsqr}. When the regressors are collinear the coefficients of the
#' dependent columns are set to zero instead of stopping.
#' 
#' @param y vector or matrix of dependent variables.
#' @param x matrix of regressors.
#' 
#' @return Matrix of coefficients, one column per column of \code{y}.
#' 
olsqr <- function(y, x){
  x <- as.matrix(x); y <- as.matrix(y)
  b <- qr.coef(qr(x), y)
  b[is.na(b)] <- 0
  return(as.matrix(b))
}

#' @title Quadratic spectral kernel
#'
#' @description The quadratic spectral kernel of Andrews (1991), used by
#' \code{jhatpr()} to weight autocovariances.
#'
#' @param x argument of the kernel (lag divided by bandwidth).
#'
#' @return Kernel weight.
#'
#' @references Andrews, D. W. K. (1991), "Heteroskedasticity and Autocorrelation Consistent Covariance Matrix Estimation", \emph{Econometrica}, 59, 817-858.
#'
#' @keywords internal
kern <- function(x){
  # procedure to evaluate the quadratic kernel at some value x.
  del <- 6*pi*x/5
  ker <- 3*(sin(del)/del-cos(del))/(del*del)
  return(ker)
}


#' @title Automatic bandwidth for the long-run covariance
#'
#' @description Andrews' (1991) automatic bandwidth for the quadratic spectral
#' kernel, from an AR(1) approximation to each column of \code{vhat} with equal
#' weights.
#'
#' @param vhat matrix whose long-run covariance is being estimated, one series per column.
#'
#' @return Bandwidth.
#'
#' @references Andrews, D. W. K. (1991), "Heteroskedasticity and Autocorrelation Consistent Covariance Matrix Estimation", \emph{Econometrica}, 59, 817-858.
#'
#' @keywords internal
bandw <- function(vhat){
  # procedure that compute the automatic bandwidth based on AR(1)
  # approximation for each vector of the matrix vhat. Each are given equal
  # weight 1.
  nt <- nrow(vhat)
  d <- ncol(vhat)
  a2n <- 0
  a2d <- 0
  for (i in 1:d){
    b <- olsqr(vhat[(2:nt),i], vhat[(1:(nt-1)),i])
    sig <- (vhat[(2:nt),i] - b%*%vhat[(1:(nt-1)),i])%*%t(vhat[(2:nt),i]-b%*%vhat[(1:(nt-1)),i])
    sig <- sig/(nt-1)
    a2n <- a2n + 4*b*b*sig*sig/(1-b)^8
    a2d <- a2d + sig*sig/(1-b)^4
  }
  a2 <- a2n/a2d
  st <- 1.3221*(a2*nt)^0.2
  return(st)
}

#' @title Long-run covariance matrix
#'
#' @description Kernel estimate of the long-run covariance matrix of the rows
#' of \code{vmat}, with the quadratic spectral kernel and the automatic
#' bandwidth of \code{bandw()}, scaled by \code{1 / (nt - d)}.
#'
#' @param vmat \code{nt} by \code{d} matrix, one series per column.
#'
#' @return \code{d} by \code{d} matrix.
#'
#' @references Andrews, D. W. K. (1991), "Heteroskedasticity and Autocorrelation Consistent Covariance Matrix Estimation", \emph{Econometrica}, 59, 817-858.
#'
#' @keywords internal
jhatpr <- function(vmat){
  # procedure to compute the long-run covariance matrix of vmat.  
  nt <- nrow(vmat)
  d <- ncol(vmat)
  jhat <- matrix(0, d, d)
  # calling the automatic bandwidth selection
  st <- bandw(vmat)
  # lag 0 covariance
  jhat <- t(vmat)%*%vmat
  # forward sum
  for (j in 1:(nt-1)){
    jhat <- jhat + c(kern(j/st))*t(matrix(vmat[((j+1):nt),], nt-j, ncol(vmat)))%*%matrix(vmat[(1:(nt-j)),], nt-j, ncol(vmat))
  }
  # backward sum
  for (j in 1:(nt-1)){
    jhat <- jhat + c(kern(j/st))*t(matrix(vmat[(1:(nt-j)),], nt-j, ncol(vmat)))%*%matrix(vmat[((j+1):nt),], nt-j, ncol(vmat))
  }
  # small sample correction
  jhat <- jhat/(nt-d)
  
  return(jhat)
}

#' @title HAC long-run covariance, with optional prewhitening
#'
#' @description Long-run covariance matrix of \code{res} by \code{jhatpr()},
#' optionally after prewhitening with a VAR(1) and recolouring (Andrews and
#' Monahan 1992). Used by the sup-F tests and the confidence intervals when
#' the errors may be serially correlated.
#'
#' @param res matrix, one series per column.
#' @param prewhit 1 to prewhiten with a VAR(1), 0 otherwise.
#'
#' @return Square matrix of the same dimension as \code{ncol(res)}.
#'
#' @references Andrews, D. W. K. (1991), "Heteroskedasticity and Autocorrelation Consistent Covariance Matrix Estimation", \emph{Econometrica}, 59, 817-858.
#' @references Andrews, D. W. K. and J. C. Monahan (1992), "An Improved Heteroskedasticity and Autocorrelation Consistent Covariance Matrix Estimator", \emph{Econometrica}, 60, 953-966.
#'
correct <- function(res,prewhit){
  d <- ncol(res)
  nt <- nrow(res)
  b <- matrix(0,d,1)
  bmat <- matrix(0,d,d)
  vstar <- matrix(0,nt-1,d)
  vmat <- res
  
  # Procedure that applies prewhitenning to the matrix vmat by filtering
  # with a VAR(1). if prewhit=0, it is skipped.
  
  if (prewhit == 1){
    #  estimating the VAR(1) 
    for (i in 1:d){
      b <- olsqr(vmat[2:nt,i],vmat[1:nt-1,])
      bmat[,i] <- b
      vstar[,i] <- vmat[2:nt,i] - vmat[1:(nt-1),]%*%b
    }
    
    # call the kernel on the residuals 
    jh <- jhatpr(vstar)
    
    # recolor 
    hac <- solve(diag(d)-bmat)%*%jh%*%t(solve(diag(d)-bmat))
  }else{
  hac <- jhatpr(vmat)
  }
  return(hac)
}




#' @title Sup-LR test against a fixed number of breaks
#'
#' @description The likelihood-ratio statistic of Qu and Perron (2007),
#' \code{2 (lr_m - lr_0)}, comparing the model with \code{m} breaks at the
#' estimated dates with the model without breaks, whose coefficients and
#' covariance are estimated by iterated feasible GLS. Serial correlation in the
#' errors is not allowed for; see \code{supft()} for that case.
#'
#' @param maty dependent variables stacked by observation: element
#' \code{(t-1)*n + i} is observation \code{t} of equation \code{i} (from \code{transf()}).
#' @param matx regressors stacked the same way, one column per coefficient (from \code{transf()}).
#' @param m number of breaks.
#' @param n number of equations.
#' @param br the estimated break dates, as an \code{m} by 1 matrix (from \code{est()}).
#' @param rbeta the stacked regime coefficients under the alternative (from \code{est()}).
#' @param bigt number of observations.
#' @param brbeta 1 if the regression coefficients change under the alternative, 0 otherwise.
#' @param brv 1 if the covariance matrix changes under the alternative, 0 otherwise.
#'
#' @return The value of the statistic. Compare with \code{critical_values()}.
#'
#' @references Qu, Z. and P. Perron (2007), "Estimating and Testing Structural Changes in Multivariate Regressions", \emph{Econometrica}, 75, 459-502.
#'
suplr <- function(maty,matx,m,n,br,rbeta,bigt,brbeta,brv){

  pmatx <- pzbar(matx,m,br,bigt)
  res <- maty - pmatx%*%rbeta
  res <- reshapeg(res,nrow(res)/n,n)
  seg <- c(0,br,(nrow(maty)/n))

  if (brv==1 && brbeta==1){
    ulr <- 0
    for (k in 1:(m+1)){
      i <- seg[k] + 1
      j <- seg[k+1]
      vvar <- t(as.matrix(res[i:j,]))%*%as.matrix(res[i:j,])/(j-i+1)
      ulr <- ulr - ((j-i+1)*n/2)*(log(2*pi)+1)-(j-i+1)/2*log(det(vvar))
    }
  }else if (brv==1 & brbeta!=1){
    # now, change in the var only, iterative method is used for the
    # optimization, with the unrestricted estimate as the prior input.
    # Then the beta is updated and then var. The procedure is iterated
    # until convergence. Non-linear optimization algorithm could be used
    # in the place of the iterative procedure
    ulr <- 0
    for (k in 1:(m+1)){
      i <- seg[k] + 1
      j <- seg[k+1]
      vvar <- t(as.matrix(res[i:j,]))%*%as.matrix(res[i:j,])/(j-i+1)
      ulr <- ulr - ((j-i+1)*n/2)*(log(2*pi)+1)-(j-i+1)/2*log(det(vvar))
    }
  }else if (brv!=1 & brbeta==1){
    ulr <- 0
    vvar <- t(as.matrix(res[1:bigt,]))%*%as.matrix(res[1:bigt,])/(bigt)
    ulr <- -(bigt*n/2)*(log(2*pi)+1)-(bigt)/2*log(det(vvar))
  }
  # now estimate the model without breaks by iterated FGLS
  b <- pracma::pinv(t(matx)%*%matx)%*%t(matx)%*%maty
  res <- maty - matx%*%b
  umat <- reshapeg(res,nrow(res)/n,n)
  vvar <- t(umat)%*%umat/bigt
  vstar <- vvar + 1
  itr <- 1
  while (max(abs(vecr(vvar-vstar)))>1e-6 & itr<1000){
    vstar <- vvar
    b <- gls_solve(matx, maty, n, pracma::pinv(vstar))
    res <- maty - matx%*%b
    umat <- reshapeg(res,nrow(res)/n,n)
    vvar <- t(umat)%*%umat/bigt
	  itr <- itr + 1
	  if (itr==1000){
	    warning(' warning: the FGLS fail, the result does not converge')
	  }
  }
  rlr <- -(bigt*n/2)*(log(2*pi)+1)-(bigt)/2*log(det(vvar))
  retp <- 2*(ulr-rlr)
  return(retp)
}


#' @title Sup-F test against a fixed number of breaks, robust to serial correlation
#'
#' @description Wald-type statistic for changes in the regression coefficients
#' at the estimated break dates, with a HAC covariance matrix from
#' \code{correct()}, so that serially correlated errors are allowed. The
#' covariance matrix of the errors is held fixed across regimes.
#'
#' @inheritParams suplr
#' @param prewhit 1 to prewhiten the scores with a VAR(1) in the HAC estimate, 0 otherwise.
#'
#' @return The value of the statistic. Compare with \code{critical_values()}.
#'
#' @references Bai, J. and P. Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes", \emph{Econometrica}, 66, 47-78.
#' @references Qu, Z. and P. Perron (2007), "Estimating and Testing Structural Changes in Multivariate Regressions", \emph{Econometrica}, 75, 459-502.
#'
supft <- function(maty,matx,m,n,br,rbeta,bigt,brbeta,brv,prewhit){
  pmatx <- pzbar(matx,m,br,bigt)
  res <- maty - pmatx%*%rbeta
  umat <- reshapeg(res, nrow(res)/n,n)
  sigma <- t(umat)%*%umat/bigt
  # (I_T %x% A) M for A = sigma^{-1/2}, block by block, without the T n by T n matrix
  isqrt <- pracma::pinv(sqrm(sigma))
  xbstar <- block_multiply(isqrt, pmatx, n)
  ubstar <- block_multiply(isqrt, res, n)
  ubstar <- reshapeg(ubstar,nrow(ubstar)/n,n)
  tempmat <- matrix(0, bigt, ncol(xbstar))
  for (k in 1:bigt){
    tempmat[k,] <- t(t(xbstar[((k-1)*n+1):((k)*n),])%*%as.matrix(ubstar[k,]))
  }
  vmat <- correct(tempmat,prewhit)
  vbeta <- pracma::pinv(t(xbstar)%*%xbstar/bigt)%*%(vmat/bigt)%*%pracma::pinv(t(xbstar)%*%xbstar/bigt)
  rsub <- matrix(0,m,m+1)
  for (j in 1:m){
    rsub[j,j] <- -1
    rsub[j,j+1] <- 1
  }
  rmat <- kronecker(rsub,diag(ncol(matx)))
  ftest <- t(rbeta)%*%t(rmat)%*%solve(rmat%*%vbeta%*%t(rmat))%*%rmat%*%rbeta
  ftest <- (bigt-(m+1)*ncol(matx))%*%ftest/(bigt)
  return(ftest)
}

#' @title Restricted sup-F test, robust to serial correlation
#'
#' @description As \code{supft()}, for a model whose coefficients satisfy the
#' linear restrictions \code{beta = R delta}.
#'
#' @inheritParams supft
#' @param R restriction matrix, \code{beta = R delta}.
#'
#' @return The value of the statistic.
#'
#' @references Bai, J. and P. Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes", \emph{Econometrica}, 66, 47-78.
#' @references Qu, Z. and P. Perron (2007), "Estimating and Testing Structural Changes in Multivariate Regressions", \emph{Econometrica}, 75, 459-502.
#'
suprft <- function(maty,matx,m,n,br,rbeta,R,bigt,brbeta,brv,prewhit){
  # local pmatx,res,xbstar,ubstar,ostar,sigma,umat,bigsig,tempmat,k,vmat,vbeta,vdelta;
  # local rmat,rsub,j,ftest;
  
  pmatx <- pzbar(matx,m,br,bigt)
  res <- maty - pmatx%*%rbeta
  umat <- reshapeg(res,nrow(res)/n,n)
  sigma <- t(umat)%*%umat/bigt
  isqrt <- pracma::pinv(sqrm(sigma))
  pmatx <- pmatx%*%R
  xbstar <- block_multiply(isqrt, pmatx, n)
  ubstar <- block_multiply(isqrt, res, n)
  ubstar <- reshapeg(ubstar,nrow(ubstar)/n,n)
  tempmat <- matrix(0,bigt,ncol(xbstar))
  for (k in 1:bigt){
    tempmat[k,] <- t(t(xbstar[((k-1)*n+1):((k)*n),])%*%as.matrix(ubstar[k,]))
  }
  vmat <- correct(tempmat,prewhit)
  vdelta <- pracma::pinv(t(xbstar)%*%xbstar/bigt)%*%(vmat/bigt)%*%pracma::pinv(t(xbstar)%*%xbstar/bigt)
  vbeta <- R%*%vdelta%*%t(R)
  rsub <- matrix(0,m,m+1)
  for (j in 1:m){
    rsub[j,j] <- -1
    rsub[j,j+1] <- 1
  }
  rmat <- kronecker(rsub, diag(ncol(matx)))
  ftest <- t(rbeta)%*%t(rmat)%*%pracma::pinv(rmat%*%vbeta%*%t(rmat))%*%rmat%*%rbeta
  ftest <- (bigt-(m+1)*ncol(matx))*ftest/(bigt)
  return(ftest)
}


#' @title GLS coefficients with a common covariance matrix
#'
#' @description Solves \code{t(X) (I \%x\% S) X b = t(X) (I \%x\% S) y} for stacked
#' data, where \code{S} is applied to every observation's block of \code{n}
#' rows, summing over observations instead of forming the \code{T n} by
#' \code{T n} matrix.
#'
#' @param X regressors stacked by observation.
#' @param y dependent variables stacked by observation.
#' @param n number of equations.
#' @param S the \code{n} by \code{n} weight matrix (an inverse covariance matrix).
#'
#' @return Coefficient vector.
#'
#' @keywords internal
gls_solve <- function(X, y, n, S){
  XS <- block_multiply(S, X, n)
  pracma::pinv(t(XS) %*% X) %*% (t(XS) %*% y)
}

#' @title Multiply every block of rows by a matrix
#'
#' @description Computes \code{(I_T \%x\% A) M} for a matrix \code{M} with
#' \code{T n} rows stacked by observation, applying the \code{n} by \code{n}
#' matrix \code{A} to each block of \code{n} rows.
#'
#' @param A \code{n} by \code{n} matrix.
#' @param M matrix with \code{T n} rows.
#' @param n number of equations.
#'
#' @return Matrix of the same shape as \code{M}.
#'
#' @keywords internal
block_multiply <- function(A, M, n){
  M <- as.matrix(M)
  bigt <- nrow(M) / n
  out <- matrix(0, nrow(M), ncol(M))
  rows <- lapply(seq_len(n), function(i) seq(i, nrow(M), by = n))
  for (i in seq_len(n)) for (k in seq_len(n)) if (A[i, k] != 0) out[rows[[i]], ] <- out[rows[[i]], ] + A[i, k] * M[rows[[k]], , drop = FALSE]
  out
}

#' @title Sequential sup-LR test, Seq(l+1 | l)
#' 
#' @description For each of the \code{l+1} segments defined by the \code{l}
#' break dates in \code{dt}, computes the sup-LR test for one additional break,
#' and returns the largest value (Section 6.3 of Qu and Perron 2007). Within a
#' segment the minimum regime length is \code{round(trm * length)}, raised to
#' the floor of \code{regime_floor()} when that is larger; segments too short
#' to hold an additional break are skipped.
#' 
#' @inheritParams est
#' @param dt the \code{l} break dates under the null, from the global search.
#' @param nseg number of segments, \code{l + 1}.
#' @param h minimum segment length used in the global search.
#' @param R_ls optional list of restriction matrices, one per segment.
#' 
#' @return The value of the test, or \code{NA} with a warning when no segment
#' is long enough to hold another break (the original code returns 0).
#' 
#' @references Qu, Z. and P. Perron (2007), "Estimating and Testing Structural Changes in Multivariate Regressions", \emph{Econometrica}, 75, 459-502.
#' 
seqlr <- function(maty,matz,n,dt,nseg,bigt,h,brbeta,brv, trm, R_ls=NULL){
  lrtest <- matrix(0,nseg,1)
  dv <- matrix(0,nseg+1,1)
  dv[2:nseg,1] <- dt
  dv[nseg+1,1] <- bigt
  floor <- regime_floor(matz, n)
  ni <- 0
  for (is in 1:nseg){
    length <- dv[is+1,1]-dv[is,1]
    hs <- max(trim_length(trm, length), floor)
    if (length >= 2*h && length >= 3*hs){
      starti <- dv[is,1]+1
      endi <- dv[is+1,1]
      segy <- as.matrix(maty[(n*(starti-1)+1):(n*endi),])
      segz <- as.matrix(matz[(n*(starti-1)+1):(n*endi),])
      if(is.null(R_ls)){
        R_tmp <- diag(ncol(matz)*2)
      }else{
        R_tmp <- R_ls[[is]]
      }
      if (brbeta==0 & brv==1){
        R_tmp <- kronecker(matrix(1,2,1),diag(ncol(matz)))
      }
      estout <- est(segy,segz,n,1,length,trm,R_tmp,brbeta,brv,h=hs)
      lrtest[is,1] <- suplr(segy,segz,1,n,estout$dx,estout$rbeta,length,brbeta,brv)
    }else{
      ni <- ni + 1
      lrtest[is,1] <- 0
    }
  }
  if (ni == nseg){
    warning("seqlr: no segment of the partition with ", nseg - 1, " break(s) is long enough to hold another break; Seq(",
            nseg, "|", nseg - 1, ") is NA.", call. = FALSE)
    return(NA_real_)
  }
  maxf <- max(lrtest[1:nseg,1])
  return(maxf)
}


#' @title Sequential sup-F test, Seq(l+1 | l), robust to serial correlation
#' 
#' @description As \code{seqlr()}, with the HAC-based sup-F test of
#' \code{supft()} in each segment. Only changes in the coefficients are
#' considered.
#' 
#' @inheritParams seqlr
#' @param prewhit 1 to prewhiten the scores with a VAR(1) before the HAC estimate, 0 otherwise.
#' 
#' @return The value of the test, or \code{NA} with a warning when no segment
#' is long enough to hold another break.
#' 
seqf <- function(maty,matz,n,dt,nseg,bigt,h,brbeta,brv,prewhit, trm, R_ls=NULL){
  lrtest <- matrix(0,nseg,1)
  dv <- matrix(0,nseg+1,1)
  dv[2:nseg,1] <- dt
  dv[nseg+1,1] <- bigt
  floor <- regime_floor(matz, n, prewhit = (prewhit == 1))
  ni <- 0
  for (is in 1:nseg){
    length <- dv[is+1,1]-dv[is,1]
    hs <- max(trim_length(trm, length), floor)
    if (length >= 2*h && length >= 3*hs){
      starti <- dv[is,1]+1
      endi <- dv[is+1,1]
      segy <- as.matrix(maty[(n*(starti-1)+1):(n*endi),])
      segz <- as.matrix(matz[(n*(starti-1)+1):(n*endi),])
      if (is.null(R_ls)){
        estout <- est(segy,segz,n,1,length,trm,diag(ncol(matz)*2),1,0,h=hs)
        lrtest[is,1] <- supft(segy,segz,1,n,estout$dx,estout$rbeta,length,1,0,prewhit)
      }else{
        estout <- est(segy,segz,n,1,length,trm,R_ls[[is]],1,0,h=hs)
        lrtest[is,1] <- suprft(segy,segz,1,n,estout$dx,estout$rbeta,R_ls[[is]],length,1,0,prewhit)
      }
    }else{
      ni <- ni + 1
      lrtest[is,1] <- 0
    }
  }
  if (ni == nseg){
    warning("seqf: no segment of the partition with ", nseg - 1, " break(s) is long enough to hold another break; Seq(",
            nseg, "|", nseg - 1, ") is NA.", call. = FALSE)
    return(NA_real_)
  }
  maxf <- max(lrtest[1:nseg,1])
  return(maxf)
}


#' @title Distribution function of the break-date estimate
#'
#' @description Cumulative distribution function of the limiting distribution
#' of the break-date estimate under shrinking shifts (Bai 1997; Theorem 3 of
#' Qu and Perron 2007), evaluated at \code{x}. The parameters are those
#' prepared by \code{cvg()}, which handles the change of scale on either side
#' of the break.
#'
#' @param x point of evaluation.
#' @param bet,alph,b,deld,gam parameters of the two branches of the distribution.
#'
#' @return Probability.
#'
#' @references Bai, J. (1997), "Estimation of a Change Point in Multiple Regression Models", \emph{Review of Economics and Statistics}, 79, 551-563.
#' @references Qu, Z. and P. Perron (2007), "Estimating and Testing Structural Changes in Multivariate Regressions", \emph{Econometrica}, 75, 459-502.
#'
#' @importFrom stats pnorm
#' @keywords internal
funcg <- function(x,bet,alph,b,deld,gam){
  if (x <= 0){
    xb <- bet*sqrt(abs(x))
    if (abs(xb) <= 30){
      g <- -sqrt(-x/(2*pi))*exp(x/8)-(bet/alph)*exp(-alph*x)*pnorm(-bet*sqrt(abs(x)))+((2*bet*bet/alph)-2-x/2)*pnorm(-sqrt(abs(x))/2)         
    }else{
      aa <- log(bet/alph)-alph*x-xb^2/2-log(sqrt(2*pi))-log(xb)
      g <- -sqrt(-x/(2*pi))*exp(x/8)-exp(aa)*pnorm(-sqrt(abs(x))/2)+((2*bet*bet/alph)-2-x/2)*pnorm(-sqrt(abs(x))/2)
    }
  }else{
    xb <- deld*sqrt(x)
    if (abs(xb) <= 30){
      g <- 1+(b/sqrt(2*pi))*sqrt(x)*exp(-b*b*x/8)+(b*deld/gam)*exp(gam*x)*pnorm(-deld*sqrt(x))+(2-b*b*x/2-2*deld*deld/gam)*pnorm(-b*sqrt(x)/2)
    }else{
      aa <- log((b*deld/gam))+gam*x-xb^2/2-log(sqrt(2*pi))-log(xb)
      g <- 1+(b/sqrt(2*pi))*sqrt(x)*exp(-b*b*x/8)+exp(aa)+(2-b*b*x/2-2*deld*deld/gam)*pnorm(-b*sqrt(x)/2) 
    }
  }
  return(g)
}



#' @title Quantiles of the break-date distribution
#'
#' @description Quantiles of the limiting distribution of a break-date
#' estimate under shrinking shifts (Bai 1997; Qu and Perron 2007, Theorem 3),
#' obtained by root finding on the distribution function \code{funcg()}. The
#' three inputs summarise the change across the break: \code{eta} is the
#' ratio of the drift parameters of the two regimes, and \code{phi1s} and
#' \code{phi2s} the ratios of long-run variance to drift on each side. The
#' search bracket widens until it contains every requested quantile; when the
#' inputs are not finite and positive the result is \code{NA}.
#'
#' @param eta ratio of the second-regime to the first-regime drift.
#' @param phi1s,phi2s scale parameters of the first and second regime.
#' @param probs probabilities of the quantiles.
#'
#' @return Vector of quantiles, one per element of \code{probs}.
#'
#' @references Bai, J. (1997), "Estimation of a Change Point in Multiple Regression Models", \emph{Review of Economics and Statistics}, 79, 551-563.
#' @references Qu, Z. and P. Perron (2007), "Estimating and Testing Structural Changes in Multivariate Regressions", \emph{Econometrica}, 75, 459-502.
#'
#' @keywords internal
cvg <- function(eta, phi1s, phi2s, probs = c(0.025, 0.05, 0.95, 0.975)){
  if (any(!is.finite(c(eta, phi1s, phi2s))) || eta <= 0 || phi1s <= 0 || phi2s <= 0){
    return(rep(NA_real_, length(probs)))
  }
  a <- phi1s/phi2s
  gam <- ((phi2s/phi1s)+1)*eta/2
  b <- sqrt(phi1s*eta/phi2s)
  deld <- sqrt(phi2s*eta/phi1s)+b/2
  alph <- a*(1+a)/2
  bet <- (1+2*a)/2
  g <- function(x) funcg(x, bet, alph, b, deld, gam)
  # widen the bracket until it contains every quantile, stopping where the
  # distribution function can no longer be evaluated
  lwb <- -200; upb <- max(200, 100/gam)
  k <- 0
  while (k < 60 && is.finite(g(2*lwb)) && g(lwb) > min(probs)){ lwb <- 2*lwb; k <- k + 1 }
  k <- 0
  while (k < 60 && is.finite(g(2*upb)) && g(upb) < max(probs)){ upb <- 2*upb; k <- k + 1 }
  gl <- g(lwb); gu <- g(upb)
  vapply(probs, function(pr){
    if (!is.finite(gl) || !is.finite(gu) || gl >= pr || gu <= pr) return(NA_real_)
    stats::uniroot(function(x) g(x) - pr, c(lwb, upb), tol = 1e-10)$root
  }, numeric(1))
}


#' @title Response-surface critical values of the sup-LR test
#'
#' @description Critical values of the sup-LR (sup-F) test against \code{k}
#' breaks from the response-surface regressions of Bai and Perron (2003),
#' valid for 1 to 10 changing parameters, on the scale of the sup-LR statistic
#' of this package. For other cases use \code{critical_values()}.
#'
#' @param k number of breaks.
#' @param q number of parameters that change at each break.
#' @param trm trimming, as a fraction of the sample.
#'
#' @return Vector of critical values at the 90, 95, 97.5 and 99 percent levels.
#'
#' @references Bai, J. and P. Perron (2003), "Critical Values for Multiple Structural Change Tests", \emph{Econometrics Journal}, 6, 72-78.
#'
cvlr <- function(k,q,trm){
  cv90 <- (7.551+1.718*q-0.041*q*q-0.610*k-15.846*trm+0.025*(q/trm))*exp(0.338/k-0.014/(trm*k))
  cv95 <- (8.238+1.756*q-0.043*q*q-0.659*k-15.436*trm+0.025*(q/trm))*exp(0.389/k-0.013/(trm*k))
  cv975 <- (8.968+1.788*q-0.045*q*q-0.715*k-15.255*trm+0.025*(q/trm))*exp(0.426/k-0.013/(trm*k))
  cv99 <- (9.879+1.771*q-0.042*q*q-0.777*k-14.551*trm+0.025*(q/trm))*exp(0.471/k-0.012/(trm*k))
  retp <- k*c(cv90,cv95,cv975,cv99)
  return(retp)
}


#' @title Response-surface critical values of the WDmax test
#'
#' @description Critical values of the WDmax test for up to five breaks from
#' the response-surface regressions of Bai and Perron (2003), valid for 1 to 10
#' changing parameters. For other cases use \code{critical_values()}.
#'
#' @inheritParams cvlr
#'
#' @return Vector of critical values at the 90, 95, 97.5 and 99 percent levels.
#'
#' @references Bai, J. and P. Perron (2003), "Critical Values for Multiple Structural Change Tests", \emph{Econometrics Journal}, 6, 72-78.
#'
cvwd <- function(q,trm){
  cv90 <- (7.316+3.128*q-8.624*trm)*exp(-0.029*q-0.412*trm)
  cv95 <- (9.039+3.318*q-9.969*trm)*exp(-0.030*q-0.327*trm)
  cv975 <- (10.703+3.465*q-11.119*trm)*exp(-0.031*q-0.250*trm)
  cv99 <- (13.189+3.346*q-11.870*trm)*exp(-0.026*q-0.173*trm)
  retp <- c(cv90,cv95,cv975,cv99)
  return(retp)
}

#' @title Response-surface critical values of the sequential test
#'
#' @description Critical values of the Seq(l+1 | l) test from the
#' response-surface regressions of Bai and Perron (2003), valid for 1 to 10
#' changing parameters. For other cases use \code{critical_values()}.
#'
#' @param l number of breaks under the null.
#' @inheritParams cvlr
#'
#' @return Vector of critical values at the 90, 95, 97.5 and 99 percent levels.
#'
#' @references Bai, J. and P. Perron (2003), "Critical Values for Multiple Structural Change Tests", \emph{Econometrics Journal}, 6, 72-78.
#'
supseq <- function(l,q,trm){
  cv90 <- (8.397+3.702*q-0.209*q*q+0.317*(l+1)-3.736/(l+1)-11.596*trm)*exp(-0.027*q+0.005*q*q)
  cv95 <- (9.879+4.086*q-0.232*q*q+0.322*(l+1)-3.687/(l+1)-11.931*trm)*exp(-0.039*q+0.006*q*q)
  cv975 <- (11.424+4.435*q-0.261*q*q+0.300*(l+1)-3.700/(l+1)-12.292*trm)*exp(-0.047*q+0.006*q*q)
  cv99 <- (13.073+4.954*q-0.317*q*q+0.285*(l+1)-3.419/(l+1)-12.452*trm)*exp(-0.058*q+0.008*q*q)
  retp <- c(cv90,cv95,cv975,cv99)
  return(retp)
}
