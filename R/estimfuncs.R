
#' @title Dating global break dates using log-likelihood 
#' 
#' @description Find the optimal partition allowing m breaks (unrestricted model)
#' 
#' @param maty vector dependent variables, with the (i-1)*n+j th element corresponding to the ith obseration from the jth equation included in the regression
#' @param matz matrix of independent variables, arranged in the same manner as maty. The number of columns corresponds to the number of basic parameters included in the regression
#' @param h    mimunum number of observations in each segment
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
dating_mle <- function(maty,matz,n,h,m,bigt){
  datevec <- matrix(0, m, m)    # up-triangular matrix contains the estimated break dates for break numbers from one to m
  optdat  <- matrix(0, bigt, m) # row index corresponds to the ending dates, column index corresponds to the number of breaks permitted before the ending date the cell contains the optimal final break date
  optmle  <- matrix(0, bigt, m) # same as above, the cell contains the MML corresponding to that break structure
  dvec    <- matrix(0, bigt, 1) # the index is the date after which we inserting the break point. The cell contains the corresponding MML
  global  <- matrix(0, m, 1)    # Global MML when i breaks are permitted
  
  
  bigvec  <- mlebigvec(maty, matz, h, n, bigt)
  
  
  # 
  # # ----- up to here
  # if (m == 1){
  #   parti_out <- parti_og(1,h,bigt-h,bigt,bigvec,bigt)
  #   datevec[1,1]  <- parti_out[1,2]
  #   glob[1,1]     <- parti_out[1,1]
  # }else{
  #   for (j1 in (2*h):bigt){
  #     parti_out <- parti_og(1,h,j1-h,j1,bigvec,bigt)           
  #     optmle[j1,1]  <- parti_out[1,1]  
  #     optdat[j1,1]  <- parti_out[1,2]  
  #   }
  #   glob[1,1]      <- optmle[bigt,1]
  #   datevec[1,1]   <- optdat[bigt,1]
  #   for (ib in 2:m){
  #     if (ib == m){
  #       jlast <- bigt
  #       for (jb in (ib*h):(jlast-h)){
  #         dvec[jb,1]  <- optmle[jb,ib-1] + bigvec[(jb+1)*bigt-jb*(jb+1)/2,1]
  #       }
  #       optmle[jlast,ib]      <- min(dvec[(ib*h):(jlast-h),1])
  #       minindcdvec           <- which.min(dvec[(ib*h):(jlast-h),1])
  #       optdat[jlast,ib]      <- (ib*h-1) + t(minindcdvec)
  #     }else{
  #       for (jlast in ((ib+1)*h):bigt){
  #         for (jb in (ib*h):(jlast-h)){
  #           dvec[jb,1]   <- optmle[jb,ib-1] + bigvec[jb*bigt-jb*(jb-1)/2+jlast-jb,1] 
  #         }
  #         optmle[jlast,ib] <- min(dvec[(ib*h):(jlast-h),1])
  #         minindcdvec  <- which.min(dvec[(ib*h):(jlast-h),1])
  #         optdat[jlast,ib] <- (ib*h-1) + t(minindcdvec)
  #       }
  #     }
  #     datevec[ib,ib] <- optdat[bigt,ib]
  #     for (i in 1:(ib-1)){
  #       xx <- ib-i
  #       datevec[xx,ib] <- optdat[datevec[xx+1,ib],xx]  
  #     }
  #     glob[ib,1] <- optmle[bigt,ib]
  #   }
  # }
  # glob <- -glob
  return(list(glob = glob,datevec = datevec,bigvec = bigvec))
}



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
