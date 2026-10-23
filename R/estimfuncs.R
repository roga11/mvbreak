
#' @title Vectorise a matrix
#'
#' @description Stacks the columns of a matrix into a vector (GAUSS \code{vec}).
#'
#' @param a matrix.
#'
#' @return Numeric vector.
#'
#' @keywords internal
vec <- function(a){
  v <- c(a)  
  return(v)
}


#' @title Regressor matrix partitioned by regime
#'
#' @description Builds the block-diagonal regressor matrix of the model with
#' \code{m} breaks: block \code{j} holds the rows of \code{matz} that belong to
#' regime \code{j}, so that \code{maty = pzbar \%*\% beta + u} with \code{beta}
#' the stacked regime coefficients.
#'
#' @param matz regressors stacked by observation (from \code{transf()}).
#' @param m number of breaks.
#' @param bb the \code{m} break dates, as an \code{m} by 1 matrix.
#' @param bigt number of observations.
#'
#' @return Matrix with \code{nrow(matz)} rows and \code{(m + 1) * ncol(matz)} columns.
#'
#' @references Bai, J. and P. Perron (2003), "Computation and Analysis of Multiple Structural Change Models", \emph{Journal of Applied Econometrics}, 18, 1-22.
#' @references Qu, Z. and P. Perron (2007), "Estimating and Testing Structural Changes in Multivariate Regressions", \emph{Econometrica}, 75, 459-502.
#'
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

#' @title Stack a system of equations for estimation
#'
#' @description Writes the system \code{y_t = (I_n \%x\% t(z_t)) S beta + u_t} in
#' stacked form: \code{maty} holds the dependent variables observation by
#' observation and \code{matz} the corresponding rows of \code{(I_n \%x\% t(z_t)) S}.
#' This is the input format of \code{est()} and of the test functions.
#'
#' @param y \code{T} by \code{n} matrix of dependent variables, one column per equation.
#' @param z \code{T} by \code{q} matrix of all regressors in the system, without duplication.
#' @param S \code{nq} by \code{p} selection matrix: column \code{k} picks the
#' regressor (within \code{I_n \%x\% t(z_t)}) that coefficient \code{k} multiplies, so a
#' zero row removes a regressor from an equation. \code{NULL} (the default) means
#' every regressor enters every equation, as in a VAR.
#'
#' @return List with \code{maty} (\code{nT} by 1) and \code{matz} (\code{nT} by \code{p}).
#'
#' @references Qu, Z. and P. Perron (2007), "Estimating and Testing Structural Changes in Multivariate Regressions", \emph{Econometrica}, 75, 459-502.
#'
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










#' @title Residuals of every regime's coefficients over the whole sample
#'
#' @description For each regime \code{j}, computes the residuals of the whole
#' sample under the coefficients of regime \code{j}. The result is the input
#' of \code{dating_m2()}, the re-partitioning step of the restricted
#' estimation (Section 4.2 of Qu and Perron 2007).
#'
#' @param maty dependent variables stacked by observation: element
#' \code{(t-1)*n + i} is observation \code{t} of equation \code{i} (from \code{transf()}).
#' @param matz regressors stacked the same way, one column per coefficient (from \code{transf()}).
#' @param nbeta stacked regime coefficients, \code{(m + 1) * q} by 1.
#' @param q number of coefficients per regime.
#' @param n number of equations.
#' @param m number of breaks.
#'
#' @return Matrix with \code{T * (m + 1)} rows and \code{n} columns; rows
#' \code{(j-1)*T + 1} to \code{j*T} hold the residuals under regime \code{j}.
#'
#' @references Qu, Z. and P. Perron (2007), "Estimating and Testing Structural Changes in Multivariate Regressions", \emph{Econometrica}, 75, 459-502.
#'
residseg <- function(maty,matz,nbeta,q,n,m){
  bigt <- nrow(maty)/n
  bigresvec <- matrix(0,bigt*(m+1),n)
  for (i in 1:(m+1)){
    bigresvec[((i-1)*bigt+1):(i*bigt),] <- reshapeg((maty - matz%*%as.matrix(nbeta[((i-1)*q+1):(i*q),])),bigt,n)
  }
  return(bigresvec)
}


#' @title Minimum length of a regime
#' 
#' @description The number of observations a regime must contain for the
#' model to be estimable: the largest number of regressors in any one equation
#' plus the number of equations, so that the residual covariance matrix is
#' nonsingular; one more when the scores are prewhitened in a HAC estimate.
#' 
#' @param matz regressor matrix stacked by observation (\code{matz} from \code{transf()}).
#' @param n number of equations.
#' @param prewhit logical; \code{TRUE} when a prewhitened HAC estimate will be computed.
#' 
#' @return Integer. The attribute \code{"p_eq"} holds the largest number of
#' regressors in one equation.
#' 
#' @keywords internal
regime_floor <- function(matz, n, prewhit = FALSE){
  p_eq <- 0L
  for (i in seq_len(n)){
    rows_i <- seq(i, nrow(matz), by = n)
    p_eq <- max(p_eq, sum(colSums(abs(matz[rows_i, , drop = FALSE])) > 0))
  }
  out <- p_eq + n + as.integer(prewhit)
  attr(out, "p_eq") <- p_eq
  out
}

#' @title Minimum segment length from the trimming
#' 
#' @description Converts the trimming fraction into a number of observations,
#' \code{round(trm * bigt)}, with R's rounding of halves to even. GAUSS rounds
#' halves up, so the two differ when \code{trm * bigt} is a half-integer; pass
#' \code{h} directly to reproduce GAUSS exactly. Every function in the package
#' obtains \code{h} from this helper.
#' 
#' @param trm trimming, as a fraction of the sample.
#' @param bigt number of observations.
#' 
#' @return Integer.
#' 
#' @keywords internal
trim_length <- function(trm, bigt){
  as.integer(round(trm * bigt))
}

#' @title Checks on the trimming and the number of breaks
#' 
#' @description Stops when the minimum segment length is below the floor
#' required by the model, or when \code{m + 1} regimes of the minimum length
#' would fill the sample so that there is no room to search; warns when a regime of the minimum length would hold fewer than
#' five observations per regressor.
#' 
#' @param h minimum segment length.
#' @param bigt number of observations.
#' @param m number of breaks.
#' @param floor value returned by \code{regime_floor()}.
#' @param fn name of the calling function, for messages.
#' @param warn logical; whether to issue the observations-per-regressor warning.
#' 
#' @return \code{NULL}, invisibly.
#' 
#' @keywords internal
check_trimming <- function(h, bigt, m, floor, fn, warn = TRUE){
  p_eq <- attr(floor, "p_eq")
  if (h < floor){
    stop(fn, ": the minimum segment length h = ", h, " is below ", as.integer(floor),
         ", the number needed to estimate a regime with ", p_eq, " regressors per equation. ",
         "Increase the trimming or the sample.", call. = FALSE)
  }
  # m + 1 regimes of at least h observations must leave room to search:
  # (m + 1) * h < T (Bai and Perron 2003: nine breaks with eps = 0.10 would fix the dates)
  m_max <- ceiling(bigt / h) - 2
  if ((m + 1) * h >= bigt){
    stop(fn, ": with T = ", bigt, " and h = ", h, " at most ", m_max,
         " break(s) can be searched; ", m, " requested.", call. = FALSE)
  }
  if (warn && p_eq > 0 && h / p_eq < 5){
    warning(fn, ": a regime of the minimum length h = ", h, " holds fewer than 5 observations per regressor (",
            round(h / p_eq, 1), "). Results in short regimes will be imprecise.", call. = FALSE)
  }
  invisible(NULL)
}

#' @title Number of workers for parallel computation
#' 
#' @description Validates the requested number of workers and caps it at the
#' number of physical cores.
#' 
#' @param workers requested number of workers; 1 means sequential.
#' @param fn name of the calling function, for messages.
#' 
#' @return Integer.
#' 
#' @keywords internal
n_workers <- function(workers, fn){
  workers <- suppressWarnings(as.integer(workers))
  if (length(workers) != 1L || is.na(workers) || workers < 1L){
    stop(fn, ": 'workers' must be a single integer of at least 1.", call. = FALSE)
  }
  cores <- parallel::detectCores(logical = FALSE)
  if (!is.na(cores) && workers > cores){
    message(fn, ": capping workers at the ", cores, " physical cores available (", workers, " requested).")
    workers <- cores
  }
  workers
}

#' @title Segment likelihoods computed in parallel
#' 
#' @description Same output as \code{mlebigvec()}: minus the log-likelihood of
#' every admissible segment, stored in the triangular layout used by
#' \code{dating_dp()}. The starting dates are split into contiguous blocks, one
#' per worker, and computed on a socket cluster that is created and stopped
#' inside the function.
#' 
#' @inheritParams dating_mle
#' @param maty dependent variables stacked by observation.
#' @param matz regressors stacked by observation.
#' @param workers number of workers; 1 computes sequentially.
#' 
#' @return Numeric vector of length \code{bigt * (bigt + 1) / 2}.
#' 
mlebigvec_par <- function(maty, matz, h, n, bigt, workers = 1){
  workers <- n_workers(workers, "mlebigvec_par")
  n_start <- bigt - h + 1
  if (workers == 1 || n_start < 2 * workers){
    return(as.numeric(mlebigvec(maty, matz, h, n, bigt)))
  }
  bounds <- round(seq(0, n_start, length.out = workers + 1))
  starts <- bounds[-length(bounds)] + 1
  ends   <- bounds[-1]
  cl <- parallel::makePSOCKcluster(workers)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  pieces <- parallel::clusterMap(cl, function(start, end, maty, matz, h, n, bigt){
    utils::getFromNamespace("mlebigvec_forpar", "mvbreak")(maty, matz, h, n, bigt, start, end)
  }, starts, ends, MoreArgs = list(maty = maty, matz = matz, h = h, n = n, bigt = bigt))
  bigvec <- unlist(pieces, use.names = FALSE)
  c(bigvec, numeric(bigt * (bigt + 1) / 2 - length(bigvec)))
}













#' @title Break dates and coefficients under restrictions
#' 
#' @description Estimates the break dates and the regime coefficients with the
#' iterative procedure of Section 4.2 of Qu and Perron (2007): the dates of the
#' unrestricted model start the iteration; the coefficients are then estimated
#' under the restrictions \code{beta = R delta}, the sample is re-partitioned
#' under those coefficients, and the two steps alternate until the dates stop
#' changing.
#' 
#' @param maty dependent variables stacked by observation (\code{maty} from \code{transf()}).
#' @param matz regressors stacked by observation (\code{matz} from \code{transf()}).
#' @param n number of equations.
#' @param m number of breaks.
#' @param bigt number of observations.
#' @param trm trimming, as a fraction of the sample. Ignored when \code{h} is given.
#' @param R restriction matrix, \code{beta = R delta}. The identity (the default) imposes no restriction.
#' @param brbeta 1 if the regression coefficients may change, 0 otherwise.
#' @param brv 1 if the covariance matrix of the errors may change, 0 otherwise.
#' @param h minimum length of a regime. By default \code{round(trm * bigt)}.
#' @param workers number of workers used to compute the segment likelihoods; 1 is sequential.
#' @param maxit largest number of re-partitioning iterations.
#' @param bigvec optional vector of segment likelihoods from \code{mlebigvec()},
#' computed once and reused when the same data are searched for several numbers
#' of breaks; \code{NULL} computes it.
#' 
#' @return List with \code{dx} (the break dates), \code{rbeta} (the coefficients
#' of all regimes, stacked) and \code{rvv} (the covariance matrices of all
#' regimes, stacked), followed by the inputs \code{maty}, \code{matz}, \code{n},
#' \code{m}, \code{bigt}, \code{h}, \code{R}, \code{brbeta} and \code{brv}, so
#' that the result can be passed to \code{confint_breaks()} and the tests.
#' 
#' @references Qu, Z. and P. Perron (2007), "Estimating and Testing Structural Changes in Multivariate Regressions", \emph{Econometrica}, 75, 459-502.
#' 
est <- function(maty, matz, n, m, bigt, trm, R = NULL, brbeta = 1, brv = 1, h = NULL,
                workers = 1, maxit = 100, bigvec = NULL){
  q <- ncol(matz)
  if (is.null(R)){
    R <- diag(q*(m+1))
  }
  floor <- regime_floor(matz, n)
  if (is.null(h)){
    h <- trim_length(trm, bigt)
    check_trimming(h, bigt, m, floor, "est")
  }else{
    check_trimming(h, bigt, m, floor, "est", warn = FALSE)
  }
  # find the optimal partition of the sample without restrictions on the
  # coefficients; this is the starting value of the iterative procedure
  if (is.null(bigvec)){
    bigvec <- mlebigvec_par(maty, matz, h, n, bigt, workers)
  } else if (length(bigvec) != bigt * (bigt + 1) / 2){
    stop("est: 'bigvec' must have T (T + 1) / 2 = ", bigt * (bigt + 1) / 2, " elements.", call. = FALSE)
  }
  datevec <- dating_dp(bigvec, h, m, bigt)$datevec
  br <- datevec[,m]
  restim_out <- r_estim(maty,matz,bigt,n,m,br,R,brbeta,brv)
  rbeta <- restim_out$beta
  rvv <- restim_out$vv
  # now find the optimal partition under the estimated coefficients, and
  # alternate until the dates stop changing
  itr <- 0
  repeat{
    bigvec2 <- residseg(maty,matz,rbeta,q,n,m)
    datevec <- dating_m2(bigvec2,h,m,n,bigt)$datevec
    if (all(datevec[,m] == br)){
      break
    }
    itr <- itr + 1
    if (itr >= maxit){
      warning("est: the break dates did not settle after ", maxit, " iterations; returning the current dates.", call. = FALSE)
      break
    }
    br <- datevec[,m]
    restim_out <- r_estim(maty,matz,bigt,n,m,br,R,brbeta,brv)
    rbeta <- restim_out$beta
    rvv <- restim_out$vv
  }
  list(dx = as.matrix(br), rbeta = rbeta, rvv = rvv,
       maty = maty, matz = matz, n = n, m = m, bigt = bigt, h = h, R = R, brbeta = brbeta, brv = brv)
}
