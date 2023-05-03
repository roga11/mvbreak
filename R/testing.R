#' @title Square root of PD matrix
#' 
#' @description compute the square root matrix of a symmetric positive definte matrix
#' 
#' @param M matrix
#' 
#' @export  
sqrm <- function(M){
  eg <- eigen(M)
  dd <- diag(eg$values)
  ve <- eg$vectors
  sq <- ve%*%sqrt(dd)%*%t(ve)
  return(sq)
}

#' @title ols parameter vector estimation
#' 
#' @details Note: This code is an adaptation of the one originally written by Yohei 
#' Yamamoto and Pierre Perron for MATLAB. Original codes can be found on 
#' Pierre Perron's website: https://blogs.bu.edu/perron/codes/
#' 
#' @references Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.
#' @references Bai, Jushan & Pierre Perron (2003), "Computation and Analysis of Multiple Structural Change Models," \emph{Journal of Applied Econometrics}, 18, 1-22.
#' 
#' @export
olsqr <- function(y, x){
  b <- solve(t(x)%*%x)%*%t(x)%*%y 
  return(as.matrix(b))
}

#' @title Compute quadratic kernel
#' 
#' @details Note: This code is an adaptation of the one originally written by Yohei 
#' Yamamoto and Pierre Perron for MATLAB. Original codes can be found on 
#' Pierre Perron's website: https://blogs.bu.edu/perron/codes/
#' 
#' @references Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.
#' @references Bai, Jushan & Pierre Perron (2003), "Computation and Analysis of Multiple Structural Change Models," \emph{Journal of Applied Econometrics}, 18, 1-22.
#' 
#' @export
kern <- function(x){
  # procedure to evaluate the quadratic kernel at some value x.
  del <- 6*pi*x/5
  ker <- 3*(sin(del)/del-cos(del))/(del*del)
  return(ker)
}


#' @title bandwidth based on AR(1) approximation
#' 
#' @details Note: This code is an adaptation of the one originally written by Yohei 
#' Yamamoto and Pierre Perron for MATLAB. Original codes can be found on 
#' Pierre Perron's website: https://blogs.bu.edu/perron/codes/
#' 
#' @references Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.
#' @references Bai, Jushan & Pierre Perron (2003), "Computation and Analysis of Multiple Structural Change Models," \emph{Journal of Applied Econometrics}, 18, 1-22.
#' 
#' @export
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

#' @title compute long run variance
#' 
#' @details Note: This code is an adaptation of the one originally written by Yohei 
#' Yamamoto and Pierre Perron for MATLAB. Original codes can be found on 
#' Pierre Perron's website: https://blogs.bu.edu/perron/codes/
#' 
#' @references Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.
#' @references Bai, Jushan & Pierre Perron (2003), "Computation and Analysis of Multiple Structural Change Models," \emph{Journal of Applied Econometrics}, 18, 1-22.
#' 
#' @export
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

#' @title HAC procedure that obtains the robust standard errors
#' 
#' @description main procedures which activates the computation of robust standard errors.
#' 
#' @references Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.
#' @references Qu & Perron (2007) Econometrica
#' 
#' @export
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




#' @title sup LRT from Qu & Perron (2007)
#' 
#' @description Procedure that constructs the supf test against m breaks. 
#' The test is the LR test under the optimal partitions computed by dating_MLE. 
#' Control variables: 
#'     if brv==1, against breaks in variance
#'     if brbeta==1, against breaks in coefficients
#'     if brv==1 and brbeta==1, allowing changes in both.
#'     
#' @note This version of the test does NOT allow for serial correlation of the errors. 
#' 
#' @references Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.
#' @references Qu & Perron (2007) Econometrica
#' 
#' @export
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
  # now estimate the restricted model
  b <- pracma::pinv(t(matx)%*%matx)%*%t(matx)%*%maty
  res <- maty - matx%*%b
  umat <- reshapeg(res,nrow(res)/n,n)
  vvar <- t(umat)%*%umat/bigt
  vstar <- vvar + 1
  itr <- 1
  while (max(abs(vecr(vvar-vstar)))>1e-6 & itr<1000){
    vstar <- vvar
    ibigv <- kronecker(diag(nrow(maty)/n), pracma::pinv(vstar))
    b <- pracma::pinv(t(matx)%*%ibigv%*%matx)%*%t(matx)%*%ibigv%*%maty 
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


#' @title supF 
#' 
#' @description Test allowing m breaks.
#' 
#' @note This version of the test DOES allow for serial correlation of the errors but does not 
#' allow for changes in the variance. 
#' 
#' @references Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.
#' @references Qu & Perron (2007) Econometrica
#' 
#' @export
supft <- function(maty,matx,m,n,br,rbeta,bigt,brbeta,brv,prewhit){
  pmatx <- pzbar(matx,m,br,bigt)
  res <- maty - pmatx%*%rbeta
  umat <- reshapeg(res, nrow(res)/n,n)
  sigma <- t(umat)%*%umat/bigt
  bigsig <- pracma::pinv(sqrm(kronecker(diag(bigt),sigma)))
  xbstar <- bigsig%*%pmatx
  ubstar <- bigsig%*%res
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

#' @title restricted supF 
#' 
#' @description restricted supF test allowing m breaks 
#' 
#' @note This version of the test DOES allow for serial correlation of the errors but does not 
#' allow for changes in the variance. 
#' 
#' @references Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.
#' @references Qu & Perron (2007) Econometrica
#' 
#' @export
suprft <- function(maty,matx,m,n,br,rbeta,R,bigt,brbeta,brv,prewhit){
  # local pmatx,res,xbstar,ubstar,ostar,sigma,umat,bigsig,tempmat,k,vmat,vbeta,vdelta;
  # local rmat,rsub,j,ftest;
  
  pmatx <- pzbar(matx,m,br,bigt)
  res <- maty - pmatx%*%rbeta
  umat <- reshapeg(res,nrow(res)/n,n)
  sigma <- t(umat)%*%umat/bigt
  bigsig <- pracma::pinv(sqrm(kronecker(diag(bigt),sigma)))
  pmatx <- pmatx%*%R
  xbstar <- bigsig%*%pmatx
  ubstar <- bigsig%*%res
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


#' @title sequential likelihood ratio 
#' 
#' @description procedure that computes the supF(l+1|l) test l=nseg-1. The l breaks
#'  under the null are taken from the global minimization (in dt).
#' 
#' @references Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.
#' @references Qu & Perron (2007) Econometrica
#' 
#' @export
seqlr <- function(maty,matz,n,dt,nseg,bigt,h,brbeta,brv, trm, R_ls=NULL){
  # remove h as input - only trm is needed but specify h as func. of trm and bigt
  # can obtain nseg from dt - consider removing as input and creating in func.
  ssr <- matrix(0,nseg,1)
  lrtest <- matrix(0,nseg,1)
  dv <- matrix(0,nseg+1,1)
  dv[2:nseg,1] <- dt
  dv[nseg+1,1] <- bigt
  ds <- matrix(0,nseg,1)
  ni <- 0
  is <- 1
  while (is <= nseg){
    length <- dv[is+1,1]-dv[is,1]
    if (length >= 2*h){
      starti <- dv[is,1]+1
      endi <- dv[is+1,1]
      segy <- as.matrix(maty[(n*(starti-1)+1):(n*endi),])
      segz <- as.matrix(matz[(n*(starti-1)+1):(n*endi),])
      if(is.null(R_ls)){
        R_tmp <- diag(ncol(matz)*2)
      }else{
        R_tmp <- R_ls[[is]]
      }
      if (brbeta==1 & brv==0){
        estout <- est(segy,segz,n,1,endi-starti+1,trm,R_tmp,1,0,h)
        lrtest[is,1] <- suplr(segy,segz,1,n,estout$dx,estout$rbeta,endi-starti+1,1,0)
      }
      if (brbeta==0 & brv==1){
        estout <- est(segy,segz,n,1,endi-starti+1,trm,kronecker(matrix(1,2,1),diag(ncol(matz))),0,1,h)
        lrtest[is,1] <- suplr(segy,segz,1,n,estout$dx,estout$rbeta,endi-starti+1,0,1)
      }
      if (brbeta==1 & brv==1){
        estout <- est(segy,segz,n,1,endi-starti+1,trm,R_tmp,1,1,h)
        lrtest[is,1] <- suplr(segy,segz,1,n,estout$dx,estout$rbeta,endi-starti+1,1,1)
      }
    }else{
      ni <- ni + 1
      lrtest[is,1] <- 0
    }
    is <- is + 1
  }
  if (ni == nseg){
    warning(paste0('Given the location of the breaks from the global optimization with ',
                   nseg-1, ' breaks there was no more place to insert an additional breaks that satisfy the minimal length requirement.'))
  }
  maxf <- max(lrtest[1:nseg,1])
  return(maxf)
}


#' @title sequential f test
#' 
#' @description procedure that computes the supF(l+1|l) test l=nseg-1. The l breaks
#'  under the null are taken from the global minimization (in dt).
#' 
#' @references Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.
#' @references Qu & Perron (2007) Econometrica
#' 
#' @export
seqf <- function(maty,matz,n,dt,nseg,bigt,h,brbeta,brv,prewhit, trm){
  ssr <- matrix(0,nseg,1)
  lrtest <- matrix(0,nseg,1)
  dv <- matrix(0,nseg+1,1)
  dv[2:nseg,1] <- dt
  dv[nseg+1,1] <- bigt
  ds <- matrix(0,nseg,1)
  ni <- 0
  is <- 1
  while (is <= nseg){
    length <- dv[is+1,1]-dv[is,1]
    if (length >= 2*h){
      starti <- dv[is,1]+1
      endi <- dv[is+1,1]
      segy <- as.matrix(maty[(n*(starti-1)+1):(n*endi),])
      segz <- as.matrix(matz[(n*(starti-1)+1):(n*endi),])
      estout <- est(segy,segz,n,1,endi-starti+1,trm,diag(ncol(matz)*2),1,0)
      lrtest[is,1] <- supft(segy,segz,1,n,estout$dx,estout$rbeta,endi-starti+1,1,0,prewhit)
    }else{
      ni <- ni + 1
      lrtest[is,1] <- 0
    }
    is <- is + 1
  }
  if (ni == nseg){
    warning(paste0('Given the location of the breaks from the global optimization with ',
                   nseg-1,' breaks there was no more place to insert an additional breaks that satisfy the minimal length requirement.')) 
  }
  maxf <- max(lrtest[1:nseg,1])
  return(maxf)
}


#' @export
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



#' @export
cvg <- function(eta,phi1s,phi2s){
  cvec <- matrix(0, 4, 1)
  a <- phi1s/phi2s
  gam <- ((phi2s/phi1s)+1)*eta/2
  b <- sqrt(phi1s*eta/phi2s)
  deld <- sqrt(phi2s*eta/phi1s)+b/2
  alph <- a*(1+a)/2
  bet <- (1+2*a)/2
  sig <- c(0.025, 0.05, 0.95, 0.975)
  isig <- 1
  while (isig <=4){
    upb <- max(c(200,(100/gam)))
    lwb <- -200
    crit <- 999999
    cct <- 1
    while (abs(crit) >= .00001){
      cct <- cct + 1
      if (cct > 1000){
        stop('the procedure to get critical values for the break dates has 
             reached the upper bound on the number of iterations. This occured
             in the procedure cvg. The resulting confidence interval for this
             break date is incorect')
      }else{
        xx <- lwb+(upb-lwb)/20
        pval <- funcg(xx,bet,alph,b,deld,gam)
        crit <- pval-sig[isig]
        if (crit <= 0){
          lwb <- xx  
        }else{
          upb <- xx
        }
      }
    }
    cvec[isig,1] <- xx
    isig <- isig + 1
  }
  return(cvec)
}


#' @title Critical values for suplr test
#' 
#' @description Critical values 
#' 
#' @references Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.
#' @references Qu & Perron (2007) Econometrica
#' @references Bai, J. and Perron, P. (2003), Critical values for multiple structural change tests. \emph{The Econometrics Journal}, 6: 72-78. https://doi.org/10.1111/1368-423X.00102
#' 
#' @export
cvlr <- function(k,q,trm){
  cv90 <- (7.551+1.718*q-0.041*q*q-0.610*k-15.846*trm+0.025*(q/trm))*exp(0.338/k-0.014/(trm*k))
  cv95 <- (8.238+1.756*q-0.043*q*q-0.659*k-15.436*trm+0.025*(q/trm))*exp(0.389/k-0.013/(trm*k))
  cv975 <- (8.968+1.788*q-0.045*q*q-0.715*k-15.255*trm+0.025*(q/trm))*exp(0.426/k-0.013/(trm*k))
  cv99 <- (9.879+1.771*q-0.042*q*q-0.777*k-14.551*trm+0.025*(q/trm))*exp(0.471/k-0.012/(trm*k))
  retp <- k*c(cv90,cv95,cv975,cv99)
  return(retp)
}


#' @title Critical values for WDmax suplr test
#' 
#' @description Critical values 
#' 
#' @references Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.
#' @references Qu & Perron (2007) Econometrica
#' @references Bai, J. and Perron, P. (2003), Critical values for multiple structural change tests. \emph{The Econometrics Journal}, 6: 72-78. https://doi.org/10.1111/1368-423X.00102
#' 
#' @export
cvwd <- function(q,trm){
  cv90 <- (7.316+3.128*q-8.624*trm)*exp(-0.029*q-0.412*trm)
  cv95 <- (9.039+3.318*q-9.969*trm)*exp(-0.030*q-0.327*trm)
  cv975 <- (10.703+3.465*q-11.119*trm)*exp(-0.031*q-0.250*trm)
  cv99 <- (13.189+3.346*q-11.870*trm)*exp(-0.026*q-0.173*trm)
  retp <- c(cv90,cv95,cv975,cv99)
  return(retp)
}

#' @title Critical values for sequential suplr test
#' 
#' @description Critical values for sequential procedure
#' 
#' @references Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.
#' @references Qu & Perron (2007) Econometrica
#' @references Bai, J. and Perron, P. (2003), Critical values for multiple structural change tests. \emph{The Econometrics Journal}, 6: 72-78. https://doi.org/10.1111/1368-423X.00102
#' 
#' @export
supseq <- function(l,q,trm){
  cv90 <- (8.397+3.702*q-0.209*q*q+0.317*(l+1)-3.736/(l+1)-11.596*trm)*exp(-0.027*q+0.005*q*q)
  cv95 <- (9.879+4.086*q-0.232*q*q+0.322*(l+1)-3.687/(l+1)-11.931*trm)*exp(-0.039*q+0.006*q*q)
  cv975 <- (11.424+4.435*q-0.261*q*q+0.300*(l+1)-3.700/(l+1)-12.292*trm)*exp(-0.047*q+0.006*q*q)
  cv99 <- (13.073+4.954*q-0.317*q*q+0.285*(l+1)-3.419/(l+1)-12.452*trm)*exp(-0.058*q+0.008*q*q)
  retp <- c(cv90,cv95,cv975,cv99)
  return(retp)
}

















