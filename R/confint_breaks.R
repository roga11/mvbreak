#' @title Confidence intervals for the break dates
#'
#' @description Confidence intervals for the estimated break dates from the
#' limiting distribution under shrinking shifts of Theorem 3 of Qu and Perron
#' (2007), which allows the coefficients, the covariance matrix of the errors,
#' or both to change. The interval for break \code{j} is
#' \code{[T_j - q_hi / a, T_j - q_lo / a]}, where \code{q_lo} and \code{q_hi} are
#' quantiles of that distribution and \code{a} is a scale built from the size
#' of the change and the second moments of the regressors and errors on each
#' side of the break, estimated from the data.
#'
#' @param fit a model fitted by \code{breaks_varp()}, \code{breaks_mvreg()},
#' \code{breaks_sur()}, \code{breaks_panel()} or \code{breaks_system()}, the
#' output of \code{est()}, or a locally ordered breaks model from
#' \code{breaks_ordered()}.
#' @param level confidence levels.
#' @param robust \code{TRUE} to allow serially correlated errors: the long-run
#' variances are then estimated with the HAC estimator of \code{correct()}.
#' @param prewhit \code{TRUE} to prewhiten with a VAR(1) in the HAC estimate.
#' Only used when \code{robust = TRUE}.
#' @param hetq \code{TRUE} to let the distribution of the regressors differ
#' across regimes, so that the second moments are estimated regime by regime.
#' With \code{hetq = FALSE} and no change in the covariance matrix they are
#' estimated from the whole sample.
#' @param B,M,seed for locally ordered breaks only: the number of replications
#' and the half-width of the grid used to simulate the distribution of
#' Theorem 4, and an optional seed (the caller's random number state is left
#' unchanged). An error is given if they are supplied for another model.
#'
#' @details For locally ordered breaks (Theorem 4) the joint distribution of
#' the two dates has no closed form: the argmax of the limiting process is
#' simulated on a grid over \code{[-M, M]} with spacing 0.02, and the bounds of
#' each date are its marginal quantiles scaled by the quadratic form \code{Pi_1j}
#' of the first group's change. The errors are assumed serially uncorrelated,
#' so \code{robust} and \code{prewhit} must be \code{FALSE}; \code{hetq} lets
#' the second moments of the regressors differ before the first break and
#' after the second. A warning is given when more than 5\% of the simulated
#' maxima fall on the edge of the grid, in which case \code{M} should be increased.
#'
#' The long-run variance of the products of standardised residuals
#' is estimated on their \code{n(n+1)/2} distinct elements and expanded with
#' the duplication matrix, which gives the same quantity as the original GAUSS
#' code when the errors are serially uncorrelated and keeps the prewhitening
#' regression well defined when they are not. Bounds are truncated to integers
#' as in the GAUSS code: the lower bound rounds down and the upper bound rounds
#' down and adds one. Unlike the GAUSS code, bounds are then clipped to the
#' sample, \code{[1, T]}: a lower bound of 1 or an upper bound of \code{T} means
#' the interval reaches the end of the sample. If the estimated scale parameters
#' are not finite and positive the bounds are \code{NA}, with a warning.
#'
#' @return A matrix with one row per break and columns \code{estimate},
#' then \code{lower} and \code{upper} for each level.
#'
#' @references Bai, J. (1997), "Estimation of a Change Point in Multiple Regression Models", \emph{Review of Economics and Statistics}, 79, 551-563.
#' @references Qu, Z. and P. Perron (2007), "Estimating and Testing Structural Changes in Multivariate Regressions", \emph{Econometrica}, 75, 459-502.
#'
#' @examples
#' y <- as.matrix(data1[, 1:2]); z <- as.matrix(data1[, 3:5])
#' fit <- breaks_sur(y, x = list(z[, 2, drop = FALSE], z[, 3, drop = FALSE]), m = 2, trim = 0.2)
#' confint_breaks(fit)
#'
#' @export
confint_breaks <- function(fit, level = c(0.90, 0.95), robust = FALSE, prewhit = FALSE, hetq = FALSE,
                           B = 10000, M = 20, seed = NULL){
  fn <- "confint_breaks"
  if (inherits(fit, "mvbreakOrdered")){
    return(confint_ordered(fit, level, robust, prewhit, hetq, B, M, seed))
  }
  if (!missing(B) || !missing(M) || !missing(seed)) stop(fn, ": 'B', 'M' and 'seed' are only used for locally ordered breaks (breaks_ordered()).", call. = FALSE)
  if (inherits(fit, "mvbreakEst")) fit <- fit$engine
  need <- c("dx", "rbeta", "rvv", "maty", "matz", "n", "m", "bigt", "brbeta", "brv")
  if (!is.list(fit) || !all(need %in% names(fit))){
    stop(fn, ": 'fit' must be a model fitted by one of the breaks_*() functions or the output of est().", call. = FALSE)
  }
  offset <- if (is.null(fit$offset)) 0L else as.integer(fit$offset)   # VAR: initial observations
  if (any(!is.finite(level)) || any(level <= 0) || any(level >= 1)){
    stop(fn, ": 'level' must lie strictly between 0 and 1.", call. = FALSE)
  }
  for (nm in c("robust", "prewhit", "hetq")){
    v <- get(nm)
    if (length(v) != 1 || is.na(v) || !is.logical(v)) stop(fn, ": '", nm, "' must be TRUE or FALSE.", call. = FALSE)
  }
  if (prewhit && !robust) warning(fn, ": 'prewhit' has no effect when robust = FALSE.", call. = FALSE)

  maty <- fit$maty; matz <- fit$matz; n <- fit$n; m <- fit$m; bigt <- fit$bigt
  p <- ncol(matz); br <- as.integer(fit$dx); beta <- as.numeric(fit$rbeta); vv <- fit$rvv
  regime_specific <- hetq || fit$brv == 1
  pw <- if (robust && prewhit) 1 else 0

  # residuals standardised by the covariance matrix of their own regime
  res <- matrix(maty - pzbar(matz, m, as.matrix(br), bigt) %*% beta, bigt, n, byrow = TRUE)
  ends <- c(0, br, bigt)
  sig <- lapply(seq_len(m + 1), function(j) vv[((j - 1) * n + 1):(j * n), , drop = FALSE])
  eta <- res
  for (j in seq_len(m + 1)){
    rows <- (ends[j] + 1):ends[j + 1]
    eta[rows, ] <- res[rows, , drop = FALSE] %*% solve(sqrm(sig[[j]]))
  }
  xblock <- function(t) matz[((t - 1) * n + 1):(t * n), , drop = FALSE]
  Dn <- duplication_matrix(n)

  out <- matrix(NA_real_, m, 1 + 2 * length(level))
  colnames(out) <- c("estimate", as.vector(rbind(paste0("lower ", 100 * level, "%"), paste0("upper ", 100 * level, "%"))))
  rownames(out) <- paste("break", seq_len(m))
  out[, 1] <- br
  for (j in seq_len(m)){
    S1 <- sig[[j]]; S2 <- sig[[j + 1]]
    iS1 <- solve(S1); iS2 <- solve(S2); S1h <- sqrm(S1); S2h <- sqrm(S2)
    A1 <- S1h %*% iS2 %*% (S2 - S1) %*% solve(S1h)
    A2 <- S2h %*% iS1 %*% (S2 - S1) %*% solve(S2h)
    delta <- beta[(j * p + 1):((j + 1) * p)] - beta[((j - 1) * p + 1):(j * p)]
    r1 <- (ends[j] + 1):ends[j + 1]; r2 <- (ends[j + 1] + 1):ends[j + 2]
    if (!regime_specific){ r1 <- seq_len(bigt); r2 <- seq_len(bigt) }

    side <- function(rows, iS_other, S_own, S_own_half){
      # second moments on one side of the break: Q, Pi and (if needed) Omega
      Q <- matrix(0, p, p); Pi <- matrix(0, p, p)
      scores <- matrix(0, length(rows), p); w <- matrix(0, length(rows), n * (n + 1) / 2)
      for (k in seq_along(rows)){
        t <- rows[k]; X <- xblock(t)
        Q <- Q + t(X) %*% iS_other %*% X
        if (!robust) Pi <- Pi + t(X) %*% iS_other %*% S_own %*% iS_other %*% X
        else scores[k, ] <- t(X) %*% iS_other %*% S_own_half %*% eta[t, ]
        if (regime_specific) w[k, ] <- vech(tcrossprod(eta[t, ]) - diag(n))
      }
      Q <- Q / length(rows)
      if (robust) Pi <- correct(scores, pw) else Pi <- Pi / length(rows)
      Om <- NULL
      if (regime_specific) Om <- if (robust) correct(w, pw) else crossprod(w) / length(rows)
      list(Q = Q, Pi = Pi, Om = Om)
    }
    s1 <- side(r1, iS2, S1, S1h); s2 <- side(r2, iS1, S2, S2h)
    Delta1 <- as.numeric(t(delta) %*% s1$Q %*% delta); Delta2 <- as.numeric(t(delta) %*% s2$Q %*% delta)
    Gam1sq <- as.numeric(t(delta) %*% s1$Pi %*% delta); Gam2sq <- as.numeric(t(delta) %*% s2$Pi %*% delta)
    if (regime_specific){
      a1 <- as.numeric(t(Dn) %*% as.vector(A1)); a2 <- as.numeric(t(Dn) %*% as.vector(A2))
      Delta1 <- Delta1 + 0.5 * sum(diag(A1 %*% A1)); Delta2 <- Delta2 + 0.5 * sum(diag(A2 %*% A2))
      Gam1sq <- Gam1sq + 0.25 * as.numeric(t(a1) %*% s1$Om %*% a1); Gam2sq <- Gam2sq + 0.25 * as.numeric(t(a2) %*% s2$Om %*% a2)
    }
    ok <- all(is.finite(c(Delta1, Delta2, Gam1sq, Gam2sq))) && Delta1 > 0 && Delta2 > 0 && Gam1sq > 0 && Gam2sq > 0
    if (!ok){
      warning(fn, ": the scale parameters for break ", j, " are not finite and positive; its bounds are NA.", call. = FALSE)
      next
    }
    probs <- as.vector(rbind((1 - level) / 2, 1 - (1 - level) / 2))
    qs <- cvg(Delta2 / Delta1, Gam1sq / Delta1, Gam2sq / Delta2, probs)
    if (any(is.na(qs))){
      warning(fn, ": the quantiles for break ", j, " could not be found; its bounds are NA.", call. = FALSE)
      next
    }
    a <- Delta1^2 / Gam1sq
    for (i in seq_along(level)){
      out[j, 2 * i]     <- max(1, trunc(br[j] - qs[2 * i] / a))
      out[j, 2 * i + 1] <- min(bigt, trunc(br[j] - qs[2 * i - 1] / a) + 1)
    }
  }
  out <- out + offset
  attr(out, "level") <- level; attr(out, "robust") <- robust; attr(out, "prewhit") <- prewhit; attr(out, "hetq") <- hetq
  out
}

#' @title Confidence intervals for locally ordered breaks
#' @description The branch of \code{confint_breaks()} for a fit of
#' \code{breaks_ordered()}: simulates the argmax of the process of Theorem 4
#' of Qu and Perron (2007) with \code{sim_ordered_argmax()} and turns its
#' marginal quantiles into bounds for the two dates.
#' @inheritParams confint_breaks
#' @return A matrix with rows \code{k1} and \code{k2}; see \code{confint_breaks()}.
#' @keywords internal
confint_ordered <- function(fit, level, robust, prewhit, hetq, B, M, seed){
  fn <- "confint_breaks"
  if (any(!is.finite(level)) || any(level <= 0) || any(level >= 1)) stop(fn, ": 'level' must lie strictly between 0 and 1.", call. = FALSE)
  for (nm in c("robust", "prewhit", "hetq")){
    v <- get(nm)
    if (length(v) != 1 || is.na(v) || !is.logical(v)) stop(fn, ": '", nm, "' must be TRUE or FALSE.", call. = FALSE)
  }
  if (robust || prewhit) stop(fn, ": the intervals for locally ordered breaks assume serially uncorrelated errors; 'robust' and 'prewhit' must be FALSE.", call. = FALSE)
  if (length(B) != 1 || is.na(B) || B < 100 || B != round(B)) stop(fn, ": 'B' must be an integer of at least 100.", call. = FALSE)
  if (length(M) != 1 || !is.finite(M) || M <= 0) stop(fn, ": 'M' must be a positive number.", call. = FALSE)
  eng <- fit[["engine"]]; n <- eng[["n"]]; bigt <- eng[["bigt"]]; p <- fit[["model"]][["p"]]
  matz <- eng[["matz"]]; k1 <- eng[["k1"]]; k2 <- eng[["k2"]]; offset <- eng[["offset"]]
  shift <- eng[["beta"]][p + seq_len(p)]; grp1 <- eng[["grp1"]]
  delt1 <- shift * grp1; delt2 <- shift * !grp1
  sigma <- crossprod(eng[["res"]]) / bigt; isig <- solve(sigma)
  Q_over <- function(rows){
    Q <- matrix(0, p, p)
    for (t in rows){ X <- matz[((t - 1) * n + 1):(t * n), , drop = FALSE]; Q <- Q + t(X) %*% isig %*% X }
    Q / length(rows)
  }
  if (hetq){
    Qj <- Q_over(seq_len(k1)); Qjp1 <- Q_over((k2 + 1):bigt)
  } else {
    Qj <- Q_over(seq_len(bigt)); Qjp1 <- Qj
  }
  form <- function(a, Q, b) as.numeric(t(a) %*% Q %*% b)
  P1j <- form(delt1, Qj, delt1); P2j <- form(delt2, Qj, delt2); P12j <- form(delt1, Qj, delt2)
  P1jp1 <- form(delt1, Qjp1, delt1); P2jp1 <- form(delt2, Qjp1, delt2); P12jp1 <- form(delt1, Qjp1, delt2)
  out <- matrix(NA_real_, 2, 1 + 2 * length(level))
  colnames(out) <- c("estimate", as.vector(rbind(paste0("lower ", 100 * level, "%"), paste0("upper ", 100 * level, "%"))))
  rownames(out) <- c("k1", "k2")
  out[, 1] <- c(k1, k2)
  if (!all(is.finite(c(P1j, P2j, P1jp1, P2jp1))) || P1j <= 0 || P2j <= 0 || P1jp1 <= 0 || P2jp1 <= 0){
    warning(fn, ": the scale parameters are not finite and positive; the bounds are NA.", call. = FALSE)
    return(out + offset)
  }
  if (!is.null(seed)){
    had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = globalenv())
    on.exit({
      if (had_seed) assign(".Random.seed", old_seed, envir = globalenv())
      else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
    }, add = TRUE)
    set.seed(seed)
  }
  v <- sim_ordered_argmax(P1j, P2j, P12j, P1jp1, P2jp1, P12jp1, M, 0.02, as.integer(B))
  edge <- mean(abs(v[, 1]) >= M - 1e-8) + mean(abs(v[, 2]) >= M - 1e-8)
  if (edge >= 0.05) warning(fn, ": ", round(100 * edge), "% of the simulated maxima lie on the edge of the grid; increase 'M'.", call. = FALSE)
  for (i in seq_along(level)){
    lo <- (1 - level[i]) / 2; hi <- 1 - lo
    for (j in 1:2){
      q_lo <- stats::quantile(v[, j], lo, names = FALSE, type = 1); q_hi <- stats::quantile(v[, j], hi, names = FALSE, type = 1)
      out[j, 2 * i]     <- max(1, trunc(out[j, 1] - q_hi / P1j))
      out[j, 2 * i + 1] <- min(bigt, trunc(out[j, 1] - q_lo / P1j) + 1)
    }
  }
  out <- out + offset
  attr(out, "level") <- level; attr(out, "robust") <- FALSE; attr(out, "prewhit") <- FALSE; attr(out, "hetq") <- hetq
  attr(out, "B") <- B; attr(out, "M") <- M
  out
}

#' @title Half-vectorisation of a symmetric matrix
#'
#' @description The distinct elements of a symmetric matrix, column by column
#' from the diagonal down.
#'
#' @param M symmetric matrix.
#'
#' @return Vector of length \code{n(n+1)/2}.
#'
#' @keywords internal
vech <- function(M){
  M[lower.tri(M, diag = TRUE)]
}

#' @title Duplication matrix
#'
#' @description The \code{n^2} by \code{n(n+1)/2} matrix \code{D} with
#' \code{vec(M) = D vech(M)} for every symmetric \code{n} by \code{n} matrix \code{M}.
#'
#' @param n dimension.
#'
#' @return Matrix of zeros and ones.
#'
#' @keywords internal
duplication_matrix <- function(n){
  D <- matrix(0, n * n, n * (n + 1) / 2)
  idx <- matrix(0, n, n); idx[lower.tri(idx, diag = TRUE)] <- seq_len(n * (n + 1) / 2)
  idx <- idx + t(idx) - diag(diag(idx))
  D[cbind(seq_len(n * n), as.vector(idx))] <- 1
  D
}
