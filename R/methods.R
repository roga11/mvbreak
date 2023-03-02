




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
  
  