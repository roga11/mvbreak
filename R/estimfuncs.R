

#' @title Model estimation
#' 
#' @description This procedure estimates the breaks and the coefficients using the restrictions, based on the iterative method
#' 
#' @param y dependent variable
#' @param z independent variable
#' @param q: total number of regression coefficients in a single regime
#' @param bigt: sample size
#' \itemize{
#'   \item{\code{}: }{. Default is \code{}.}
#'   \item{\code{}: }{. Default is \code{}.}
#' }
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
est <- function(maty,matz,n,m,bigt,trm,R,brbeta,brv){
  q <- ncol(matz)
  h <- round(trm*bigt)
  tstar <- bigt
  datemle_out <- dating_mle(maty,matz,n,h,m,bigt)
  #global <- datemle_out$global # not used 
  datevec <- datemle_out$datevec
  #bigvec <- datemle_out$bigvec # not used
  
  # find the optimal partiton of  the sample without restrictions 
  # on the coefficients, serves as the starting value of the 
  # iterative procedure@
  br <- datevec[,m] #break dates
  restim_out <- r_estim(maty,matz,bigt,n,m,datevec[,m],R,brbeta,brv)
  rbeta <- restim_out$nbeta
  #rvv <- restim_out$nvv  # not used 
  
  # now find optimal partitions under the estimated coefficients
  diff <- 0
  while (diff!=-1){
    bigvec2 <- residuals(maty,matz,rbeta,q,n,m)
    datem2_out <- dating_m2(bigvec2,h,m,n,bigt)
    #global <- datem2_out$global  # not used 
    datevec <- datem2_out$datevec
    #bigvec <- datem2_out$bigvec  # not used 
    if (datevec[,m]==br){
      diff <- -1
      dx <- br 
    }else{
      br <- datevec[,m]
      restim_out <- r_estim(maty,matz,bigt,n,m,br,R,brbeta,brv)
      rbeta <- restim_out$nbeta
      #rvv <- restim_out$nvv  # not used 
      diff <-diff + 1               
    }
  }
  return(list(dx = dx, rbeta = rbeta))
}
