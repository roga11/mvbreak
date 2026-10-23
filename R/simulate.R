# Data-generating processes with structural breaks, one function per model type.
# Each returns the data together with a description of the process (class
# "mvbreakSim") that mc_breaks() uses to redraw it and to know the true values.

#' @title Simulate a VAR with structural breaks
#'
#' @description Draws \code{y_t = c_j + A_1j y_{t-1} + ... + A_pj y_{t-p} + u_t},
#' \code{u_t ~ N(0, Sigma_j)}, with regime \code{j} in force up to the
#' \code{j}-th break date. The coefficients, the covariance matrix, or both may
#' differ across regimes.
#'
#' @param n_obs number of observations.
#' @param coef the coefficients: an \code{n} by \code{1 + n p} matrix
#' \code{[c, A_1, ..., A_p]} (one row per equation) for a process without
#' breaks, or a list of \code{m + 1} such matrices, one per regime.
#' @param sigma the covariance matrix of the errors, or a list of \code{m + 1} matrices.
#' @param dates the \code{m} true break dates (the last observation of each
#' regime but the last), or \code{NULL} for no break. A single matrix in
#' \code{coef} or \code{sigma} is used in every regime.
#' @param burnin number of initial observations discarded.
#' @param seed optional seed; the caller's random number state is left unchanged.
#'
#' @return An object of class \code{"mvbreakSim"}: a list with \code{y}
#' (\code{n_obs} by \code{n}), \code{x} (\code{NULL} here), and \code{dgp}, the
#' description of the process: \code{type}, \code{n}, \code{m}, \code{dates},
#' \code{coef} and \code{sigma} as lists over regimes, \code{coef_matrix} (the
#' true coefficients in the layout of \code{coef()} on the fitted model),
#' \code{change}, \code{lags}, and \code{args} (the arguments, so that the
#' process can be redrawn).
#'
#' @examples
#' A1 <- rbind(c(0.5, 0.1), c(0.0, 0.3)); A2 <- rbind(c(-0.2, 0.0), c(0.2, 0.6))
#' sim <- simulate_varp(200, coef = list(cbind(0, A1), cbind(1, A2)), sigma = diag(2),
#'                      dates = 100, seed = 1)
#' fit <- breaks_varp(sim$y, p = 1, m = 1)
#' fit$dates
#'
#' @export
simulate_varp <- function(n_obs, coef, sigma, dates = NULL, burnin = 100, seed = NULL){
  fn <- "simulate_varp"
  args <- list(n_obs = n_obs, coef = coef, sigma = sigma, dates = dates, burnin = burnin)
  d <- regimes_from(coef, sigma, dates, n_obs, fn)
  n <- nrow(d$coef[[1]]); k <- ncol(d$coef[[1]])
  if ((k - 1) %% n != 0 || k < 1 + n) stop(fn, ": each coefficient matrix must be n by 1 + n p.", call. = FALSE)
  p <- (k - 1) %/% n
  reg <- c(rep(1L, burnin), regime_index(d$dates, n_obs))
  u <- draw_errors(n_obs + burnin, d$sigma, reg, seed)
  y <- matrix(0, n_obs + burnin, n)
  for (t in seq_len(n_obs + burnin)){
    lagged <- if (p > 0) as.vector(t(y[pmax(t - seq_len(p), 1), , drop = FALSE][seq_len(p), , drop = FALSE])) else numeric(0)
    if (t <= p) lagged <- rep(0, n * p)
    y[t, ] <- d$coef[[reg[t]]] %*% c(1, lagged) + u[t, ]
  }
  y <- y[burnin + seq_len(n_obs), , drop = FALSE]; colnames(y) <- paste0("y", seq_len(n))
  eq <- colnames(y); zn <- c("const", paste0(rep(eq, p), ".l", rep(seq_len(p), each = n)))
  new_sim(y, NULL, type = "VAR", d = d, coef_matrix = coef_matrix_from(d$coef, eq, zn), lags = p, args = args, fn = fn)
}

#' @title Simulate a multivariate regression with structural breaks
#'
#' @description Draws \code{y_t = B_j x_t + u_t}, \code{u_t ~ N(0, Sigma_j)},
#' with the same regressors in every equation. The regressors are supplied or
#' drawn as independent standard normals.
#'
#' @inheritParams simulate_varp
#' @param coef an \code{n} by \code{k} matrix \code{B} (one row per equation,
#' the intercept first when \code{intercept = TRUE}), or a list of \code{m + 1} such matrices.
#' @param x \code{n_obs} by \code{k_x} matrix of regressors, or \code{NULL} to
#' draw \code{n_x} standard normal regressors.
#' @param n_x number of regressors to draw when \code{x = NULL}.
#' @param intercept whether the first column of \code{coef} multiplies a constant.
#'
#' @return An object of class \code{"mvbreakSim"}; see \code{simulate_varp()}.
#' \code{x} holds the regressors used.
#'
#' @examples
#' sim <- simulate_mvreg(150, coef = list(rbind(c(0, 1), c(0, -1)), rbind(c(1, 1), c(1, -1))),
#'                       sigma = diag(2), dates = 75, n_x = 1, seed = 2)
#' breaks_mvreg(sim$y, sim$x, m = 1)$dates
#'
#' @export
simulate_mvreg <- function(n_obs, coef, sigma, dates = NULL, x = NULL, n_x = 1, intercept = TRUE, seed = NULL){
  fn <- "simulate_mvreg"
  args <- list(n_obs = n_obs, coef = coef, sigma = sigma, dates = dates, x = x, n_x = n_x, intercept = intercept)
  d <- regimes_from(coef, sigma, dates, n_obs, fn)
  n <- nrow(d$coef[[1]]); k <- ncol(d$coef[[1]])
  reg <- regime_index(d$dates, n_obs)
  rs <- rng_scope(seed); on.exit(rs(), add = TRUE)
  if (is.null(x)){
    if (n_x > 0) x <- matrix(stats::rnorm(n_obs * n_x), n_obs, n_x, dimnames = list(NULL, paste0("x", seq_len(n_x)))) else x <- matrix(0, n_obs, 0)
  } else {
    x <- check_matrix(x, fn, "x"); if (nrow(x) != n_obs) stop(fn, ": 'x' must have n_obs rows.", call. = FALSE)
    if (is.null(colnames(x))) colnames(x) <- paste0("x", seq_len(ncol(x)))
  }
  if (k != intercept + ncol(x)) stop(fn, ": each coefficient matrix must have intercept + ncol(x) = ", intercept + ncol(x), " columns.", call. = FALSE)
  u <- draw_errors(n_obs, d$sigma, reg, NULL)
  X <- cbind(if (intercept) 1, x)
  y <- matrix(0, n_obs, n)
  for (t in seq_len(n_obs)) y[t, ] <- d$coef[[reg[t]]] %*% X[t, ] + u[t, ]
  colnames(y) <- paste0("y", seq_len(n))
  zn <- c(if (intercept) "const", colnames(x))
  new_sim(y, if (ncol(x)) x else NULL, type = "multivariate regression", d = d,
          coef_matrix = coef_matrix_from(d$coef, colnames(y), zn), lags = 0L, args = args, fn = fn, intercept = intercept)
}

#' @title Simulate seemingly unrelated regressions with structural breaks
#'
#' @description Draws \code{y_it = t(b_ij) x_it + u_it}, \code{u_t ~ N(0, Sigma_j)},
#' with equation-specific regressors, supplied or drawn as standard normals.
#'
#' @inheritParams simulate_mvreg
#' @param coef a list with one vector per equation (the intercept first when
#' \code{intercept = TRUE}) for a process without breaks, or a list of
#' \code{m + 1} such lists, one per regime.
#' @param x list of \code{n} regressor matrices with \code{n_obs} rows
#' (\code{NULL} entries for equations with only an intercept), or \code{NULL}
#' to draw \code{n_x[i]} standard normal regressors for equation \code{i}.
#' @param n_x number of regressors to draw for each equation when \code{x = NULL}.
#'
#' @return An object of class \code{"mvbreakSim"}; see \code{simulate_varp()}.
#' \code{x} is the list of regressors used.
#'
#' @examples
#' sim <- simulate_sur(150, coef = list(list(c(0, 1), c(0, 2)), list(c(1, 1), c(1, 2))),
#'                     sigma = diag(2), dates = 75, n_x = c(1, 1), seed = 3)
#' breaks_sur(sim$y, sim$x, m = 1)$dates
#'
#' @export
simulate_sur <- function(n_obs, coef, sigma, dates = NULL, x = NULL, n_x = 1, intercept = TRUE, seed = NULL){
  fn <- "simulate_sur"
  args <- list(n_obs = n_obs, coef = coef, sigma = sigma, dates = dates, x = x, n_x = n_x, intercept = intercept)
  d <- regimes_from(coef, sigma, dates, n_obs, fn, is_sur = TRUE)
  n <- length(d$coef[[1]])
  reg <- regime_index(d$dates, n_obs)
  rs <- rng_scope(seed); on.exit(rs(), add = TRUE)
  n_x <- rep_len(n_x, n)
  if (is.null(x)){
    x <- lapply(seq_len(n), function(i) if (n_x[i] > 0) matrix(stats::rnorm(n_obs * n_x[i]), n_obs, n_x[i], dimnames = list(NULL, paste0("x", seq_len(n_x[i])))) else NULL)
  } else {
    if (!is.list(x) || length(x) != n) stop(fn, ": 'x' must be a list with one element per equation (", n, ").", call. = FALSE)
    x <- lapply(seq_len(n), function(i){
      if (is.null(x[[i]])) return(NULL)
      xi <- check_matrix(x[[i]], fn, paste0("x[[", i, "]]")); if (nrow(xi) != n_obs) stop(fn, ": every element of 'x' must have n_obs rows.", call. = FALSE)
      if (is.null(colnames(xi))) colnames(xi) <- paste0("x", seq_len(ncol(xi))); xi })
  }
  for (i in seq_len(n)) if (length(d$coef[[1]][[i]]) != intercept + NCOL(x[[i]]) * !is.null(x[[i]]))
    stop(fn, ": the coefficient vector of equation ", i, " must have intercept + number of regressors = ", intercept + NCOL(x[[i]]) * !is.null(x[[i]]), " elements.", call. = FALSE)
  u <- draw_errors(n_obs, d$sigma, reg, NULL)
  y <- matrix(0, n_obs, n)
  for (i in seq_len(n)){
    Xi <- cbind(if (intercept) rep(1, n_obs), x[[i]])
    for (t in seq_len(n_obs)) y[t, i] <- sum(Xi[t, ] * d$coef[[reg[t]]][[i]]) + u[t, i]
  }
  colnames(y) <- paste0("y", seq_len(n))
  coef_names <- unlist(lapply(seq_len(n), function(i) paste(colnames(y)[i], c(if (intercept) "const", colnames(x[[i]])), sep = ":")))
  cm <- sapply(d$coef, function(cj) unlist(cj)); dimnames(cm) <- list(coef_names, paste("regime", seq_along(d$coef)))
  new_sim(y, x, type = "SUR", d = d, coef_matrix = cm, lags = 0L, args = args, fn = fn, intercept = intercept)
}

#' @title Simulate a panel with common structural breaks
#'
#' @description Draws \code{y_it = a_ij + t(b_j) x_t + u_it} for units
#' \code{i = 1..N}: unit-specific intercepts, slopes common to all units,
#' regressors common to all units (supplied or drawn), errors \code{N(0, Sigma_j)}
#' correlated across units, and breaks common to every unit.
#'
#' @inheritParams simulate_mvreg
#' @param intercepts vector of \code{N} unit intercepts, or a list of \code{m + 1} such vectors.
#' @param slopes vector of \code{k} pooled slopes, or a list of \code{m + 1} such vectors.
#'
#' @return An object of class \code{"mvbreakSim"}; see \code{simulate_varp()}.
#'
#' @examples
#' sim <- simulate_panel(150, intercepts = list(c(0, 0, 0), c(1, 1, 1)), slopes = list(0.5, -0.5),
#'                       sigma = diag(3), dates = 75, n_x = 1, seed = 4)
#' breaks_panel(sim$y, sim$x, m = 1)$dates
#'
#' @export
simulate_panel <- function(n_obs, intercepts, slopes, sigma, dates = NULL, x = NULL, n_x = 1, seed = NULL){
  fn <- "simulate_panel"
  args <- list(n_obs = n_obs, intercepts = intercepts, slopes = slopes, sigma = sigma, dates = dates, x = x, n_x = n_x)
  m <- length(dates)
  as_list <- function(v, what){ if (is.list(v)){ if (length(v) != m + 1) stop(fn, ": '", what, "' must have m + 1 = ", m + 1, " elements.", call. = FALSE); v } else rep(list(v), m + 1) }
  intercepts <- as_list(intercepts, "intercepts"); slopes <- as_list(slopes, "slopes")
  N <- length(intercepts[[1]]); k <- length(slopes[[1]])
  coef <- lapply(seq_len(m + 1), function(j) cbind(intercepts[[j]], matrix(slopes[[j]], N, k, byrow = TRUE)))
  out <- simulate_mvreg(n_obs, coef = coef, sigma = sigma, dates = dates, x = x, n_x = n_x, intercept = TRUE, seed = seed)
  units <- paste0("unit", seq_len(N)); colnames(out$y) <- units
  cm <- sapply(seq_len(m + 1), function(j) c(intercepts[[j]], slopes[[j]]))
  dimnames(cm) <- list(c(paste(units, "const", sep = ":"), colnames(out$x)), paste("regime", seq_len(m + 1)))
  out$dgp$type <- "panel"; out$dgp$coef_matrix <- cm; out$dgp$args <- args; out$dgp$fn <- fn
  out
}

#' @title Simulate locally ordered breaks
#'
#' @description Draws a system with equation-specific regressors in which the
#' coefficients of the equations in \code{first} change after \code{k1} and those
#' of the other equations after \code{k2 >= k1}, with one covariance matrix
#' (the model of Section 5 of Qu and Perron 2007 and of \code{breaks_ordered()}).
#'
#' @inheritParams simulate_sur
#' @param coef_before,coef_after lists with one coefficient vector per equation
#' (the intercept first when \code{intercept = TRUE}), before and after each
#' equation's own break.
#' @param sigma covariance matrix of the errors.
#' @param k1,k2 the break dates of the first and second group, \code{k1 <= k2}.
#' @param first the equations (indices) whose break comes first.
#'
#' @return An object of class \code{"mvbreakSim"} with \code{type = "ordered"}
#' and \code{dates = c(k1, k2)}; see \code{simulate_varp()}.
#'
#' @examples
#' sim <- simulate_ordered(120, coef_before = list(c(0, 1), c(0, 1)),
#'                         coef_after = list(c(1, 1), c(1, 1)), sigma = diag(2),
#'                         k1 = 50, k2 = 56, n_x = c(1, 1), seed = 5)
#' breaks_ordered(sim$y, sim$x)$dates
#'
#' @export
simulate_ordered <- function(n_obs, coef_before, coef_after, sigma, k1, k2, first = 1, x = NULL, n_x = 1, intercept = TRUE, seed = NULL){
  fn <- "simulate_ordered"
  args <- list(n_obs = n_obs, coef_before = coef_before, coef_after = coef_after, sigma = sigma, k1 = k1, k2 = k2, first = first, x = x, n_x = n_x, intercept = intercept)
  n <- length(coef_before)
  if (length(coef_after) != n) stop(fn, ": 'coef_before' and 'coef_after' must have one vector per equation.", call. = FALSE)
  if (k1 > k2 || k1 < 1 || k2 >= n_obs) stop(fn, ": need 1 <= k1 <= k2 < n_obs.", call. = FALSE)
  first <- sort(unique(as.integer(first)))
  if (any(first < 1) || any(first > n) || length(first) == n) stop(fn, ": 'first' must index a strict subset of the equations.", call. = FALSE)
  rs <- rng_scope(seed); on.exit(rs(), add = TRUE)
  n_x <- rep_len(n_x, n)
  if (is.null(x)){
    x <- lapply(seq_len(n), function(i) if (n_x[i] > 0) matrix(stats::rnorm(n_obs * n_x[i]), n_obs, n_x[i], dimnames = list(NULL, paste0("x", seq_len(n_x[i])))) else NULL)
  } else {
    if (!is.list(x) || length(x) != n) stop(fn, ": 'x' must be a list with one element per equation (", n, ").", call. = FALSE)
    x <- lapply(seq_len(n), function(i){
      if (is.null(x[[i]])) return(NULL)
      xi <- check_matrix(x[[i]], fn, paste0("x[[", i, "]]")); if (nrow(xi) != n_obs) stop(fn, ": every element of 'x' must have n_obs rows.", call. = FALSE)
      if (is.null(colnames(xi))) colnames(xi) <- paste0("x", seq_len(ncol(xi))); xi })
  }
  u <- draw_errors(n_obs, list(sigma), rep(1L, n_obs), NULL)
  y <- matrix(0, n_obs, n)
  for (i in seq_len(n)){
    Xi <- cbind(if (intercept) rep(1, n_obs), x[[i]])
    if (ncol(Xi) != length(coef_before[[i]]) || ncol(Xi) != length(coef_after[[i]])) stop(fn, ": the coefficient vectors of equation ", i, " must have ", ncol(Xi), " elements.", call. = FALSE)
    kb <- if (i %in% first) k1 else k2
    for (t in seq_len(n_obs)) y[t, i] <- sum(Xi[t, ] * (if (t <= kb) coef_before[[i]] else coef_after[[i]])) + u[t, i]
  }
  colnames(y) <- paste0("y", seq_len(n))
  coef_names <- unlist(lapply(seq_len(n), function(i) paste(colnames(y)[i], c(if (intercept) "const", colnames(x[[i]])), sep = ":")))
  cm <- cbind(before = unlist(coef_before), after = unlist(coef_after)); rownames(cm) <- coef_names
  d <- list(coef = list(coef_before, coef_after), sigma = list(sigma), dates = c(k1 = k1, k2 = k2), m = 1L)
  out <- new_sim(y, x, type = "ordered", d = d, coef_matrix = cm, lags = 0L, args = args, fn = fn, intercept = intercept)
  out$dgp$first <- first; out$dgp$change <- "coef"
  out
}

# ---- helpers ------------------------------------------------------------------

#' @title Regime lists from the coefficient and covariance arguments
#' @description Expands single matrices to one per regime, checks lengths against
#' the break dates, and decides what changes across regimes.
#' @param coef,sigma,dates,n_obs the arguments of the simulation functions.
#' @param fn name of the calling function.
#' @param is_sur whether \code{coef} is a list of per-equation vectors (SUR layout).
#' @return List with \code{coef}, \code{sigma} (lists over regimes), \code{dates}, \code{m}, \code{change}.
#' @keywords internal
regimes_from <- function(coef, sigma, dates, n_obs, fn, is_sur = FALSE){
  if (!is.null(dates)){
    dates <- as.integer(dates)
    if (any(diff(dates) <= 0) || dates[1] < 1 || dates[length(dates)] >= n_obs) stop(fn, ": 'dates' must be increasing and inside 1..n_obs - 1.", call. = FALSE)
  }
  m <- length(dates)
  per_regime <- if (is_sur) is.list(coef) && is.list(coef[[1]]) else is.list(coef)
  coef <- if (per_regime) coef else rep(list(coef), m + 1)
  sigma <- if (is.list(sigma)) sigma else rep(list(sigma), m + 1)
  if (length(coef) != m + 1) stop(fn, ": 'coef' must have m + 1 = ", m + 1, " elements, one per regime.", call. = FALSE)
  if (length(sigma) != m + 1) stop(fn, ": 'sigma' must have m + 1 = ", m + 1, " elements, one per regime.", call. = FALSE)
  coef_change <- m > 0 && any(vapply(coef[-1], function(cj) !isTRUE(all.equal(cj, coef[[1]])), logical(1)))
  cov_change <- m > 0 && any(vapply(sigma[-1], function(sj) !isTRUE(all.equal(sj, sigma[[1]])), logical(1)))
  change <- if (coef_change && cov_change) "both" else if (cov_change) "cov" else "coef"
  list(coef = coef, sigma = sigma, dates = dates, m = m, change = change)
}

#' @title Draw normal errors regime by regime
#' @keywords internal
#' @param n_obs number of draws.
#' @param sigma list of covariance matrices.
#' @param reg regime index of each observation.
#' @param seed optional seed, scoped to this call.
#' @return \code{n_obs} by \code{n} matrix.
draw_errors <- function(n_obs, sigma, reg, seed){
  rs <- rng_scope(seed); on.exit(rs(), add = TRUE)
  n <- nrow(sigma[[1]])
  e <- matrix(stats::rnorm(n_obs * n), n_obs, n)
  u <- matrix(0, n_obs, n)
  for (j in seq_along(sigma)){
    rows <- which(reg == j)
    if (length(rows)) u[rows, ] <- e[rows, , drop = FALSE] %*% chol(sigma[[j]])
  }
  u
}

#' @title Scoped seeding
#' @description Sets the seed if one is given and returns a function that
#' restores the caller's random number state.
#' @param seed optional seed.
#' @return A function to call on exit.
#' @keywords internal
rng_scope <- function(seed){
  if (is.null(seed)) return(function() NULL)
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv()) else NULL
  set.seed(seed)
  function(){
    if (had) assign(".Random.seed", old, envir = globalenv())
    else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
  }
}

#' @title True coefficients in the layout of the fitted model
#' @keywords internal
#' @param coef list of \code{n} by \code{k} matrices, one per regime.
#' @param eq equation names.
#' @param zn regressor names.
#' @return \code{n k} by \code{m + 1} matrix with the names \code{coef()} uses.
coef_matrix_from <- function(coef, eq, zn){
  cm <- sapply(coef, function(B) as.vector(t(B)))
  cm <- matrix(cm, ncol = length(coef))
  dimnames(cm) <- list(as.vector(t(outer(eq, zn, paste, sep = ":"))), paste("regime", seq_along(coef)))
  cm
}

#' @title Build a simulated-data object
#' @keywords internal
#' @param y,x the data.
#' @param type model type.
#' @param d output of \code{regimes_from()} (or the equivalent for ordered breaks).
#' @param coef_matrix true coefficients in the fitted layout.
#' @param lags number of lags (VAR).
#' @param args the generator's arguments, without the seed.
#' @param fn the generator's name.
#' @param intercept whether the model has an intercept.
#' @return An object of class \code{"mvbreakSim"}.
new_sim <- function(y, x, type, d, coef_matrix, lags, args, fn, intercept = TRUE){
  structure(list(y = y, x = x,
                 dgp = list(type = type, n = ncol(y), m = d$m, dates = d$dates, coef = d$coef, sigma = d$sigma,
                            coef_matrix = coef_matrix, change = d$change, lags = lags, intercept = intercept,
                            n_obs = nrow(y), args = args, fn = fn)),
            class = "mvbreakSim")
}

#' @title Print a simulated-data object
#' @param x an object of class \code{"mvbreakSim"}.
#' @param ... unused.
#' @return \code{x}, invisibly.
#' @export
print.mvbreakSim <- function(x, ...){
  d <- x[["dgp"]]
  cat("Simulated ", d[["type"]], if (d[["lags"]] > 0) paste0("(", d[["lags"]], ")"), " with ", d[["n"]], " equations and ", d[["n_obs"]], " observations\n", sep = "")
  if (d[["type"]] == "ordered") cat("Locally ordered breaks at k1 = ", d[["dates"]][["k1"]], ", k2 = ", d[["dates"]][["k2"]], "\n", sep = "")
  else if (d[["m"]] == 0) cat("No break\n")
  else cat(d[["m"]], " break(s) at ", paste(d[["dates"]], collapse = ", "), "; change in the ",
           switch(d[["change"]], coef = "coefficients", cov = "covariance matrix", both = "coefficients and covariance matrix"), "\n", sep = "")
  invisible(x)
}
