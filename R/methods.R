
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
  In  <- diag(ncol(y))
  maty <- as.matrix(y[1,])
  if (is.null(S)==FALSE){
    matz <- kronecker(In, t(z[1,]))%*%S 
  }else{
    matz <- kronecker(In, t(z[1,]))
  }
  for (i in 2:nrow(y)){
    maty <- rbind(maty, as.matrix(y[i,]))
    
    if (is.null(S)==FALSE){
      matz_tmp <- kronecker(In, t(z[i,]))%*%S 
    }else{
      matz_tmp <- kronecker(In, t(z[i,]))
    }
    matz <- rbind(matz, matz_tmp)
  } 
  return(list(maty = maty, matz = matz))
}
  
  