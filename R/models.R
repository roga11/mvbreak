# User-facing model functions. Each one builds the stacked system (y, z, S) and the
# cross-regime restriction matrix R for one model type, hands them to the engine
# est(), and returns a classed object that confint_breaks() and the tests take.

#' @title Structural changes in a system of equations
#'
#' @description The general entry point: the model
#' \code{y_t = (I_n \%x\% t(z_t)) S beta_j + u_t} with \code{m} breaks, where the
#' selection matrix \code{S} says which regressors enter which equation and the
#' restriction matrix \code{R} which coefficients are held constant across
#' regimes. The model-specific functions \code{breaks_varp()},
#' \code{breaks_mvreg()}, \code{breaks_sur()} and \code{breaks_panel()} build
#' \code{z}, \code{S} and \code{R} from plain arguments and call this.
#'
#' @param y \code{T} by \code{n} matrix of dependent variables, one column per equation.
#' @param z \code{T} by \code{q} matrix of all regressors in the system, without duplication.
#' @param S \code{nq} by \code{p} selection matrix (see \code{transf()}); \code{NULL}
#' means every regressor enters every equation.
#' @param m number of breaks to estimate.
#' @param trim minimum length of a regime as a fraction of the sample.
#' @param change what changes at the breaks: the regression coefficients
#' (\code{"coef"}), the covariance matrix of the errors (\code{"cov"}) or both.
#' @param fixed coefficients held constant across regimes, as indices into the
#' \code{p} coefficients or as names (see \code{coef_names} in the result);
#' \code{NULL} lets every coefficient change. Ignored when \code{change = "cov"}.
#' @param R restriction matrix \code{beta = R delta}, built from \code{fixed} when
#' \code{NULL}. Supplying it overrides \code{fixed}.
#' @param h minimum length of a regime in observations; by default \code{round(trim * T)}.
#' @param workers number of workers for the segment likelihoods; 1 is sequential.
#' @param coef_names optional names for the \code{p} coefficients.
#' @param equation_names optional names for the equations; by default the column names of \code{y}.
#'
#' @return An object of class \code{"mvbreakEst"}: a list with
#' \item{dates}{the estimated break dates;}
#' \item{coef}{a \code{p} by \code{m + 1} matrix of coefficients, one column per regime;}
#' \item{sigma}{a list of the \code{m + 1} error covariance matrices;}
#' \item{model}{the model type, the number of equations and the names of equations and coefficients;}
#' \item{settings}{\code{m}, \code{trim}, \code{h}, \code{change}, \code{fixed}, \code{S}, \code{R} and
#' \code{R_user} (whether \code{R} was supplied rather than built from \code{fixed});}
#' \item{data}{\code{y} and \code{z};}
#' \item{engine}{the output of \code{est()}, used by \code{confint_breaks()} and the tests;}
#' \item{call}{the call.}
#'
#' @references Qu, Z. and P. Perron (2007), "Estimating and Testing Structural Changes in Multivariate Regressions", \emph{Econometrica}, 75, 459-502.
#'
#' @examples
#' # the two-equation SUR of the Qu and Perron example: equation 1 uses regressors 1 and 2,
#' # equation 2 uses regressors 1 and 3; two breaks in the coefficients
#' y <- as.matrix(data1[, 1:2]); z <- as.matrix(data1[, 3:5])
#' S <- rbind(c(1,0,0,0), c(0,1,0,0), c(0,0,0,0), c(0,0,1,0), c(0,0,0,0), c(0,0,0,1))
#' fit <- breaks_system(y, z, S, m = 2, trim = 0.2, change = "coef")
#' fit$dates
#' # only the intercepts (coefficients 1 and 3) change
#' fit_r <- breaks_system(y, z, S, m = 2, trim = 0.2, change = "coef", fixed = c(2, 4))
#' fit_r$dates
#'
#' @export
breaks_system <- function(y, z, S = NULL, m = 1, trim = 0.15, change = c("coef", "cov", "both"),
                          fixed = NULL, R = NULL, h = NULL, workers = 1,
                          coef_names = NULL, equation_names = NULL){
  fn <- if (is.null(sys.call(-1)) || !grepl("^breaks_", deparse(sys.call(-1)[[1]]))) "breaks_system" else deparse(sys.call(-1)[[1]])
  change <- match.arg(change)
  y <- check_matrix(y, fn, "y"); z <- check_matrix(z, fn, "z")
  if (nrow(z) != nrow(y)) stop(fn, ": 'y' and 'z' must have the same number of rows.", call. = FALSE)
  n <- ncol(y); q <- ncol(z); bigt <- nrow(y)
  if (!is.null(S)){
    S <- check_matrix(S, fn, "S")
    if (nrow(S) != n * q) stop(fn, ": 'S' must have n * q = ", n * q, " rows.", call. = FALSE)
    if (any(colSums(abs(S)) == 0)) stop(fn, ": every column of 'S' must select at least one regressor.", call. = FALSE)
  }
  p <- if (is.null(S)) n * q else ncol(S)
  if (length(m) != 1 || is.na(m) || m < 1 || m != round(m)) stop(fn, ": 'm' must be a positive integer.", call. = FALSE)
  if (length(trim) != 1 || is.na(trim) || trim <= 0 || trim >= 0.5) stop(fn, ": 'trim' must lie strictly between 0 and 0.5.", call. = FALSE)
  if (is.null(equation_names)) equation_names <- if (!is.null(colnames(y))) colnames(y) else paste0("y", seq_len(n))
  if (is.null(coef_names)) coef_names <- default_coef_names(n, q, S, equation_names, colnames(z))
  if (length(coef_names) != p) stop(fn, ": 'coef_names' must have length p = ", p, ".", call. = FALSE)

  brbeta <- as.integer(change != "cov"); brv <- as.integer(change != "coef")
  R_user <- !is.null(R) && change != "cov"
  if (change == "cov"){
    if (!is.null(fixed) || !is.null(R)) message(fn, ": 'fixed' and 'R' are not used when only the covariance matrix changes.")
    R <- kronecker(matrix(1, m + 1, 1), diag(p)); fixed_idx <- integer(0)
  } else {
    fixed_idx <- fixed_indices(fixed, coef_names, fn)
    if (is.null(R)){
      R <- restriction_matrix(p, m, setdiff(seq_len(p), fixed_idx))
    } else {
      R <- check_matrix(R, fn, "R")
      if (nrow(R) != p * (m + 1)) stop(fn, ": 'R' must have p * (m + 1) = ", p * (m + 1), " rows.", call. = FALSE)
    }
  }
  stacked <- transf(y, z, S)
  # the trimming checks are run here so that messages name the user-facing function
  if (is.null(h)) h <- trim_length(trim, bigt)
  check_trimming(h, bigt, m, regime_floor(stacked$matz, n), fn)
  engine <- est(stacked$maty, stacked$matz, n, m, bigt, trim, R, brbeta, brv, h = h, workers = workers)
  coef <- matrix(engine$rbeta, p, m + 1, dimnames = list(coef_names, paste("regime", seq_len(m + 1))))
  sigma <- lapply(seq_len(m + 1), function(j){
    sg <- engine$rvv[((j - 1) * n + 1):(j * n), , drop = FALSE]; dimnames(sg) <- list(equation_names, equation_names); sg })
  names(sigma) <- paste("regime", seq_len(m + 1))
  structure(list(dates = as.integer(engine$dx), coef = coef, sigma = sigma,
                 model = list(type = "system", n = n, equations = equation_names, coef_names = coef_names, p = p),
                 settings = list(m = m, trim = trim, h = engine$h, change = change, fixed = fixed_idx, S = S, R = R, R_user = R_user),
                 data = list(y = y, z = z), engine = engine, call = match.call()),
            class = "mvbreakEst")
}

#' @title Structural changes in a VAR
#'
#' @description A vector autoregression of order \code{p} with intercept,
#' \code{y_t = c_j + A_1j y_{t-1} + ... + A_pj y_{t-p} + u_t}, with \code{m} common
#' breaks in the coefficients, the error covariance matrix, or both. The first
#' \code{p} observations are used as initial values.
#'
#' @inheritParams breaks_system
#' @param p lag order.
#' @param fixed coefficients held constant across regimes, by index or by name;
#' names are \code{"<equation>:const"} and \code{"<equation>:<variable>.l<lag>"}.
#' @param fixed_equations equations (by index or name) whose coefficients are all held constant.
#'
#' @return An object of class \code{"mvbreakEst"}; see \code{breaks_system()}.
#'
#' @examples
#' y <- as.matrix(data1[, 1:2]); colnames(y) <- c("y1", "y2")
#' fit <- breaks_varp(y, p = 1, m = 2, trim = 0.2, change = "coef")
#' fit$dates
#' fit$coef
#'
#' @export
breaks_varp <- function(y, p = 1, m = 1, trim = 0.15, change = c("coef", "cov", "both"),
                        fixed = NULL, fixed_equations = NULL, h = NULL, workers = 1){
  fn <- "breaks_varp"
  y <- check_matrix(y, fn, "y")
  if (length(p) != 1 || is.na(p) || p < 1 || p != round(p)) stop(fn, ": 'p' must be a positive integer.", call. = FALSE)
  bigt <- nrow(y); n <- ncol(y)
  if (bigt <= p + 1) stop(fn, ": too few observations for ", p, " lag(s).", call. = FALSE)
  eq <- if (!is.null(colnames(y))) colnames(y) else paste0("y", seq_len(n))
  lags <- do.call(cbind, lapply(seq_len(p), function(l) y[(p + 1 - l):(bigt - l), , drop = FALSE]))
  z <- cbind(1, lags); zn <- c("const", paste0(rep(eq, p), ".l", rep(seq_len(p), each = n))); colnames(z) <- zn
  yy <- y[(p + 1):bigt, , drop = FALSE]; colnames(yy) <- eq
  coef_names <- as.vector(t(outer(eq, zn, paste, sep = ":")))
  fixed <- c(fixed_indices(fixed, coef_names, fn), equation_indices(fixed_equations, eq, coef_names, fn))
  out <- breaks_system(yy, z, S = NULL, m = m, trim = trim, change = change, fixed = if (length(fixed)) fixed else NULL,
                       h = h, workers = workers, coef_names = coef_names, equation_names = eq)
  out$model$type <- "VAR"; out$model$lags <- p; out$dates <- out$dates + p; out$call <- match.call()
  out$engine$offset <- p
  out
}

#' @title Structural changes in a multivariate regression
#'
#' @description A system of \code{n} equations with the same regressors in every
#' equation, \code{y_t = B_j x_t + u_t}, with \code{m} common breaks. With
#' \code{x = NULL} the model is a vector of means and a covariance matrix,
#' which is the case of breaks in the mean and variance of a multivariate series.
#'
#' @inheritParams breaks_system
#' @param x \code{T} by \code{k} matrix of regressors common to all equations, or \code{NULL}.
#' @param intercept whether to include an intercept in every equation.
#' @param fixed coefficients held constant across regimes, by index or by name;
#' names are \code{"<equation>:<regressor>"}.
#' @param fixed_equations equations (by index or name) whose coefficients are all held constant.
#'
#' @return An object of class \code{"mvbreakEst"}; see \code{breaks_system()}.
#'
#' @examples
#' y <- as.matrix(data1[, 1:2]); x <- as.matrix(data1[, 4:5])
#' fit <- breaks_mvreg(y, x, m = 2, trim = 0.2, change = "both")
#' fit$dates
#'
#' @export
breaks_mvreg <- function(y, x = NULL, intercept = TRUE, m = 1, trim = 0.15, change = c("coef", "cov", "both"),
                         fixed = NULL, fixed_equations = NULL, h = NULL, workers = 1){
  fn <- "breaks_mvreg"
  y <- check_matrix(y, fn, "y"); n <- ncol(y); bigt <- nrow(y)
  eq <- if (!is.null(colnames(y))) colnames(y) else paste0("y", seq_len(n))
  if (!is.null(x)){
    x <- check_matrix(x, fn, "x")
    if (nrow(x) != bigt) stop(fn, ": 'x' must have as many rows as 'y'.", call. = FALSE)
    xn <- if (!is.null(colnames(x))) colnames(x) else paste0("x", seq_len(ncol(x)))
  } else { xn <- character(0) }
  if (!intercept && is.null(x)) stop(fn, ": the model needs an intercept or regressors.", call. = FALSE)
  z <- cbind(if (intercept) matrix(1, bigt, 1), x); zn <- c(if (intercept) "const", xn); colnames(z) <- zn
  coef_names <- as.vector(t(outer(eq, zn, paste, sep = ":")))
  fixed <- c(fixed_indices(fixed, coef_names, fn), equation_indices(fixed_equations, eq, coef_names, fn))
  out <- breaks_system(y, z, S = NULL, m = m, trim = trim, change = change, fixed = if (length(fixed)) fixed else NULL,
                       h = h, workers = workers, coef_names = coef_names, equation_names = eq)
  out$model$type <- "multivariate regression"; out$call <- match.call()
  out
}

#' @title Structural changes in seemingly unrelated regressions
#'
#' @description A system of \code{n} equations with equation-specific regressors,
#' \code{y_it = t(b_ij) x_it + u_it}, with \code{m} common breaks.
#'
#' @inheritParams breaks_system
#' @param x list of \code{n} matrices of regressors, one per equation, each with
#' \code{T} rows, or \code{NULL} for an equation with only an intercept.
#' @param intercept whether to include an intercept in every equation.
#' @param fixed coefficients held constant across regimes, by index or by name;
#' names are \code{"<equation>:<regressor>"}.
#' @param fixed_equations equations (by index or name) whose coefficients are all held constant.
#'
#' @return An object of class \code{"mvbreakEst"}; see \code{breaks_system()}.
#'
#' @examples
#' # the Qu and Perron example: equation 1 uses regressors 1 and 2, equation 2 uses 1 and 3,
#' # where regressor 1 is a constant; two breaks in the coefficients
#' y <- as.matrix(data1[, 1:2]); z <- as.matrix(data1[, 3:5])
#' fit <- breaks_sur(y, x = list(z[, 2, drop = FALSE], z[, 3, drop = FALSE]), m = 2, trim = 0.2)
#' fit$dates
#' # only the intercepts change
#' fit_r <- breaks_sur(y, x = list(z[, 2, drop = FALSE], z[, 3, drop = FALSE]), m = 2, trim = 0.2,
#'                     fixed = c("y1:x1", "y2:x1"))
#' fit_r$dates
#'
#' @export
breaks_sur <- function(y, x, intercept = TRUE, m = 1, trim = 0.15, change = c("coef", "cov", "both"),
                       fixed = NULL, fixed_equations = NULL, h = NULL, workers = 1){
  fn <- "breaks_sur"
  y <- check_matrix(y, fn, "y")
  d <- build_sur(y, x, intercept, fn)
  fixed <- c(fixed_indices(fixed, d$coef_names, fn), equation_indices(fixed_equations, d$eq, d$coef_names, fn))
  out <- breaks_system(y, d$z, d$S, m = m, trim = trim, change = change, fixed = if (length(fixed)) fixed else NULL,
                       h = h, workers = workers, coef_names = d$coef_names, equation_names = d$eq)
  out$model$type <- "SUR"; out$call <- match.call()
  out
}

#' @title Build the stacked regressors and selection matrix of a SUR
#' @description Shared by \code{breaks_sur()} and \code{breaks_ordered()}: puts
#' the intercept and every equation's regressors in one matrix \code{z} without
#' duplication and builds the selection matrix \code{S} and the coefficient names.
#' @param y \code{T} by \code{n} matrix of dependent variables (already checked).
#' @param x list of \code{n} regressor matrices, \code{NULL} entries allowed.
#' @param intercept whether every equation has an intercept.
#' @param fn name of the calling function, for messages.
#' @return List with \code{z}, \code{S}, \code{coef_names}, \code{eq} (equation
#' names) and \code{equation} (the equation of each coefficient).
#' @keywords internal
build_sur <- function(y, x, intercept, fn){
  n <- ncol(y); bigt <- nrow(y)
  eq <- if (!is.null(colnames(y))) colnames(y) else paste0("y", seq_len(n))
  if (!is.list(x) || length(x) != n) stop(fn, ": 'x' must be a list with one element per equation (", n, ").", call. = FALSE)
  z <- if (intercept) matrix(1, bigt, 1) else matrix(0, bigt, 0); zn <- if (intercept) "const" else character(0)
  sel <- vector("list", n); coef_names <- character(0); equation <- integer(0)
  for (i in seq_len(n)){
    cols <- if (intercept) 1L else integer(0)
    if (!is.null(x[[i]])){
      xi <- check_matrix(x[[i]], fn, paste0("x[[", i, "]]"))
      if (nrow(xi) != bigt) stop(fn, ": every element of 'x' must have as many rows as 'y'.", call. = FALSE)
      xin <- if (!is.null(colnames(xi))) colnames(xi) else paste0("x", seq_len(ncol(xi)))
      cols <- c(cols, ncol(z) + seq_len(ncol(xi))); z <- cbind(z, xi); zn <- c(zn, xin)
    }
    if (length(cols) == 0) stop(fn, ": equation ", i, " has no regressors.", call. = FALSE)
    sel[[i]] <- cols; coef_names <- c(coef_names, paste(eq[i], zn[cols], sep = ":")); equation <- c(equation, rep(i, length(cols)))
  }
  colnames(z) <- zn; q <- ncol(z); p <- length(coef_names)
  S <- matrix(0, n * q, p); k <- 0
  for (i in seq_len(n)) for (cc in sel[[i]]){ k <- k + 1; S[(i - 1) * q + cc, k] <- 1 }
  list(z = z, S = S, coef_names = coef_names, eq = eq, equation = equation)
}

#' @title Locally ordered breaks in a system of equations
#'
#' @description The model of Section 5 of Qu and Perron (2007): two groups of
#' equations, each with one break in all of its coefficients, the first group's
#' break at \code{k1} and the second group's at \code{k2 >= k1}, close to each
#' other, with a covariance matrix of the errors that does not change. The
#' dates are estimated jointly by a grid search over every admissible pair,
#' with FGLS at each pair. The LR test of no break against the ordered breaks
#' uses the pairs with \code{k2 - k1 <= max_gap}, so that its limiting
#' distribution is that of the sup-LR test for one break with all \code{p}
#' coefficients changing; \code{test_breaks()} reports it with its critical
#' values and \code{confint_breaks()} gives the intervals of Theorem 4.
#'
#' @param y \code{T} by \code{n} matrix of dependent variables, one column per equation.
#' @param x regressors: a \code{T} by \code{k} matrix common to all equations
#' (a multivariate regression), a list of \code{n} matrices, one per equation
#' (a SUR; \code{NULL} entries for equations with only an intercept), or
#' \code{NULL}.
#' @param intercept whether to include an intercept in every equation.
#' @param lags number of lags of \code{y} added as regressors to every
#' equation, which gives a VAR when \code{x = NULL}. The first \code{lags}
#' observations are used as initial values.
#' @param first the equations (by index or name) whose break comes first; the
#' other equations break at the second date.
#' @param trim minimum length of the first and last regimes as a fraction of the sample.
#' @param max_gap largest distance \code{k2 - k1} allowed when computing the
#' test; by default \code{trunc(sqrt(T))}. The estimates are not restricted by it.
#' @param h minimum length of the first and last regimes in observations; by default \code{round(trim * T)}.
#'
#' @return An object of class \code{"mvbreakOrdered"}: a list with
#' \item{dates}{the two estimated dates, \code{k1} for the first group and \code{k2} for the second;}
#' \item{coef}{a \code{p} by 2 matrix: each coefficient before and after the break of its own equation;}
#' \item{sigma}{the covariance matrix of the errors;}
#' \item{lr}{the LR statistic of no break against ordered breaks, and the pair at which it is attained;}
#' \item{model}{the model type, the number of equations, names, and which coefficients belong to the first group;}
#' \item{settings}{\code{trim}, \code{h}, \code{max_gap}, \code{first}, \code{lags};}
#' \item{data}{\code{y} and \code{z};}
#' \item{engine}{the stacked data, the full likelihood grid and the fit at the estimated pair, used by \code{test_breaks()} and \code{confint_breaks()};}
#' \item{call}{the call.}
#'
#' @references Qu, Z. and P. Perron (2007), "Estimating and Testing Structural Changes in Multivariate Regressions", \emph{Econometrica}, 75, 459-502.
#'
#' @examples
#' # the example of the original code: equation 1 uses regressor 1, equation 2 uses regressor 2,
#' # both with an intercept; the break in equation 1 comes first
#' y <- as.matrix(data2[, 1:2]); z <- as.matrix(data2[, 3:5])
#' fit <- breaks_ordered(y, x = list(z[, 2, drop = FALSE], z[, 3, drop = FALSE]), trim = 0.2)
#' fit$dates
#' fit$lr
#'
#' @export
breaks_ordered <- function(y, x = NULL, intercept = TRUE, lags = 0, first = 1, trim = 0.15, max_gap = NULL, h = NULL){
  fn <- "breaks_ordered"
  y <- check_matrix(y, fn, "y"); n <- ncol(y)
  if (n < 2) stop(fn, ": the model needs at least two equations.", call. = FALSE)
  if (length(lags) != 1 || is.na(lags) || lags < 0 || lags != round(lags)) stop(fn, ": 'lags' must be a non-negative integer.", call. = FALSE)
  if (length(trim) != 1 || is.na(trim) || trim <= 0 || trim >= 0.5) stop(fn, ": 'trim' must lie strictly between 0 and 0.5.", call. = FALSE)
  eq <- if (!is.null(colnames(y))) colnames(y) else paste0("y", seq_len(n))
  if (is.character(first)){
    bad <- setdiff(first, eq)
    if (length(bad)) stop(fn, ": unknown equation name(s) in 'first': ", paste(bad, collapse = ", "), ".", call. = FALSE)
    first_idx <- match(first, eq)
  } else {
    first_idx <- as.integer(first)
    if (anyNA(first_idx) || any(first_idx < 1) || any(first_idx > n)) stop(fn, ": 'first' must give equation indices between 1 and ", n, " or names.", call. = FALSE)
  }
  first_idx <- sort(unique(first_idx))
  if (length(first_idx) == n) stop(fn, ": 'first' names every equation; at least one must break second.", call. = FALSE)
  # regressors per equation: common matrix, list, or none; lags of y added to every equation
  bigt0 <- nrow(y)
  if (is.null(x)) x_list <- vector("list", n)
  else if (is.list(x)){
    if (length(x) != n) stop(fn, ": 'x' must be a matrix, a list with one element per equation (", n, "), or NULL.", call. = FALSE)
    x_list <- x
  } else {
    x <- check_matrix(x, fn, "x")
    if (nrow(x) != bigt0) stop(fn, ": 'x' must have as many rows as 'y'.", call. = FALSE)
    if (is.null(colnames(x))) colnames(x) <- paste0("x", seq_len(ncol(x)))
    x_list <- replicate(n, x, simplify = FALSE)
  }
  design <- if (is.list(x)) "SUR" else if (!is.null(x)) "multivariate regression" else if (lags > 0) "VAR" else "means"
  if (lags > 0){
    if (bigt0 <= lags + 1) stop(fn, ": too few observations for ", lags, " lag(s).", call. = FALSE)
    L <- do.call(cbind, lapply(seq_len(lags), function(l) y[(lags + 1 - l):(bigt0 - l), , drop = FALSE]))
    colnames(L) <- paste0(rep(eq, lags), ".l", rep(seq_len(lags), each = n))
    x_list <- lapply(x_list, function(xi){
      if (is.null(xi)) return(L)
      xi <- check_matrix(xi, fn, "x"); if (nrow(xi) != bigt0) stop(fn, ": every element of 'x' must have as many rows as 'y'.", call. = FALSE)
      cbind(xi[(lags + 1):bigt0, , drop = FALSE], L)
    })
    y <- y[(lags + 1):bigt0, , drop = FALSE]
  }
  colnames(y) <- eq
  if (!intercept && all(vapply(x_list, is.null, logical(1)))) stop(fn, ": the model needs an intercept or regressors.", call. = FALSE)
  d <- build_sur(y, x_list, intercept, fn)
  bigt <- nrow(y); p <- ncol(d$S)
  grp1 <- d$equation %in% first_idx
  if (is.null(max_gap)) max_gap <- as.integer(trunc(sqrt(bigt)))
  if (length(max_gap) != 1 || is.na(max_gap) || max_gap < 0 || max_gap != round(max_gap)) stop(fn, ": 'max_gap' must be a non-negative integer.", call. = FALSE)
  stacked <- transf(y, d$z, d$S)
  if (is.null(h)) h <- trim_length(trim, bigt)
  check_trimming(h, bigt, 1, regime_floor(stacked$matz, n), fn)
  g <- ordered_grid(stacked$maty, stacked$matz, n, bigt, grp1, h, as.integer(max_gap))
  best <- g$best; bg <- g$best_gap
  if (is.null(bg)) stop(fn, ": no admissible pair satisfies 'max_gap'.", call. = FALSE)
  # the no-break model, for the LR statistic
  null <- fgls_null(stacked$maty, stacked$matz, n, bigt)
  ll0 <- null$loglik
  lr <- 2 * (bg$loglik - ll0)
  beta <- as.numeric(best$beta)
  coef <- cbind(beta[seq_len(p)], beta[seq_len(p)] + beta[p + seq_len(p)])
  dimnames(coef) <- list(d$coef_names, c("before", "after"))
  sigma <- best$vvar; dimnames(sigma) <- list(eq, eq)
  engine <- list(maty = stacked$maty, matz = stacked$matz, n = n, bigt = bigt, grp1 = grp1, h = h,
                 k1 = best$k1, k2 = best$k2, beta = beta, res = best$res, loglik = g$loglik,
                 k1_test = bg$k1, k2_test = bg$k2, loglik_test = bg$loglik, loglik_null = ll0, offset = as.integer(lags))
  structure(list(dates = c(k1 = best$k1 + lags, k2 = best$k2 + lags), coef = coef, sigma = sigma,
                 lr = list(statistic = lr, dates = c(k1 = bg$k1 + lags, k2 = bg$k2 + lags), max_gap = max_gap),
                 model = list(type = "ordered", design = design, n = n, equations = eq, coef_names = d$coef_names, p = p,
                              first = eq[first_idx], first_coef = grp1, lags = lags),
                 settings = list(trim = trim, h = h, max_gap = max_gap, first = first_idx, lags = lags),
                 data = list(y = y, z = d$z), engine = engine, call = match.call()),
            class = "mvbreakOrdered")
}

#' @title Structural changes in a panel with common breaks
#'
#' @description \code{N} units observed over \code{T} periods, each with its own
#' intercept and with slope coefficients common to all units (pooled),
#' \code{y_it = a_ij + t(b_j) x_it + u_it}, with \code{m} breaks common to every
#' unit. The errors may be correlated across units. The pooling is imposed
#' through the selection matrix; \code{pooled = FALSE} gives every unit its
#' own slopes, which is the SUR of \code{breaks_sur()}.
#'
#' @inheritParams breaks_system
#' @param x regressors: a \code{T} by \code{k} matrix common to all units, or a
#' list of \code{N} such matrices, one per unit (same \code{k} for every unit).
#' @param pooled whether the slopes are common to all units.
#' @param fixed coefficients held constant across regimes, by index or by name;
#' names are \code{"<unit>:const"} for the intercepts and \code{"<regressor>"}
#' for the pooled slopes.
#'
#' @return An object of class \code{"mvbreakEst"}; see \code{breaks_system()}.
#'
#' @examples
#' y <- as.matrix(data1[, 1:2]); x <- as.matrix(data1[, 4, drop = FALSE])
#' fit <- breaks_panel(y, x, m = 1, trim = 0.2)
#' fit$dates
#'
#' @export
breaks_panel <- function(y, x = NULL, pooled = TRUE, m = 1, trim = 0.15, change = c("coef", "cov", "both"),
                         fixed = NULL, h = NULL, workers = 1){
  fn <- "breaks_panel"
  y <- check_matrix(y, fn, "y"); N <- ncol(y); bigt <- nrow(y)
  units <- if (!is.null(colnames(y))) colnames(y) else paste0("unit", seq_len(N))
  if (is.null(x)){
    out <- breaks_mvreg(y, NULL, intercept = TRUE, m = m, trim = trim, change = change, fixed = fixed, h = h, workers = workers)
    out$model$type <- "panel"; out$call <- match.call(); return(out)
  }
  if (!is.list(x)) x <- replicate(N, x, simplify = FALSE)
  if (length(x) != N) stop(fn, ": 'x' must be one matrix or a list with one matrix per unit.", call. = FALSE)
  x <- lapply(seq_len(N), function(i) check_matrix(x[[i]], fn, paste0("x[[", i, "]]")))
  k <- ncol(x[[1]])
  if (any(vapply(x, ncol, 1L) != k) || any(vapply(x, nrow, 1L) != bigt)) stop(fn, ": every unit needs the same ", k, " regressors over ", bigt, " periods.", call. = FALSE)
  xn <- if (!is.null(colnames(x[[1]]))) colnames(x[[1]]) else paste0("x", seq_len(k))
  if (!pooled){
    out <- breaks_sur(y, x, intercept = TRUE, m = m, trim = trim, change = change, fixed = fixed, h = h, workers = workers)
    out$model$type <- "panel"; out$call <- match.call(); return(out)
  }
  # z holds the intercept and every unit's regressors; S gives unit i its own intercept
  # and maps its regressors to the k pooled slope coefficients
  z <- cbind(1, do.call(cbind, x)); q <- ncol(z); p <- N + k
  S <- matrix(0, N * q, p)
  for (i in seq_len(N)){
    S[(i - 1) * q + 1, i] <- 1
    for (cc in seq_len(k)) S[(i - 1) * q + 1 + (i - 1) * k + cc, N + cc] <- 1
  }
  coef_names <- c(paste(units, "const", sep = ":"), xn)
  out <- breaks_system(y, z, S, m = m, trim = trim, change = change, fixed = fixed, h = h, workers = workers,
                       coef_names = coef_names, equation_names = units)
  out$model$type <- "panel"; out$model$pooled <- TRUE; out$call <- match.call()
  out
}

# ---- helpers ------------------------------------------------------------------

#' @title Coerce to a numeric matrix with checks
#' @description Accepts a numeric vector, matrix or data frame; stops on
#' missing or non-numeric values.
#' @param x object to coerce.
#' @param fn,arg names used in messages.
#' @return Numeric matrix.
#' @keywords internal
check_matrix <- function(x, fn, arg){
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.null(dim(x))) x <- matrix(x, ncol = 1)
  if (!is.numeric(x)) stop(fn, ": '", arg, "' must be numeric.", call. = FALSE)
  if (anyNA(x)) stop(fn, ": '", arg, "' has missing values.", call. = FALSE)
  x
}

#' @title Default coefficient names
#' @description Names of the \code{p} coefficients of a stacked system, built
#' from the equation and regressor names.
#' @param n number of equations.
#' @param q number of regressors.
#' @param S selection matrix or \code{NULL}.
#' @param eq equation names.
#' @param zn regressor names, possibly \code{NULL}.
#' @return Character vector of length \code{p}.
#' @keywords internal
default_coef_names <- function(n, q, S, eq, zn){
  if (is.null(zn)) zn <- paste0("z", seq_len(q))
  full <- as.vector(t(outer(eq, zn, paste, sep = ":")))   # order of I_n %x% t(z_t): equation-major
  if (is.null(S)) return(full)
  apply(S, 2, function(col){ rows <- which(col != 0); if (length(rows) == 1) full[rows] else paste(full[rows], collapse = "=") })
}

#' @title Indices of coefficients held constant
#' @description Turns the \code{fixed} argument (indices or names) into
#' coefficient indices, with checks.
#' @param fixed indices, names or \code{NULL}.
#' @param coef_names names of the coefficients.
#' @param fn name used in messages.
#' @return Integer vector.
#' @keywords internal
fixed_indices <- function(fixed, coef_names, fn){
  if (is.null(fixed)) return(integer(0))
  p <- length(coef_names)
  if (is.character(fixed)){
    idx <- match(fixed, coef_names)
    if (anyNA(idx)) stop(fn, ": unknown coefficient name(s) in 'fixed': ", paste(fixed[is.na(idx)], collapse = ", "),
                         ". Available: ", paste(coef_names, collapse = ", "), ".", call. = FALSE)
  } else if (is.logical(fixed) && length(fixed) == p){
    idx <- which(fixed)
  } else {
    idx <- as.integer(fixed)
    if (anyNA(idx) || any(idx < 1) || any(idx > p)) stop(fn, ": 'fixed' must give coefficient indices between 1 and ", p, " or names.", call. = FALSE)
  }
  if (length(idx) == p) stop(fn, ": every coefficient is held constant; nothing can change.", call. = FALSE)
  sort(unique(idx))
}

#' @title Indices of the coefficients of whole equations
#' @description Turns \code{fixed_equations} (indices or names of equations)
#' into the indices of all their coefficients, assuming coefficient names of
#' the form \code{"<equation>:<regressor>"}.
#' @param fixed_equations indices, names or \code{NULL}.
#' @param eq equation names.
#' @param coef_names coefficient names.
#' @param fn name used in messages.
#' @return Integer vector.
#' @keywords internal
equation_indices <- function(fixed_equations, eq, coef_names, fn){
  if (is.null(fixed_equations)) return(integer(0))
  if (is.character(fixed_equations)){
    bad <- setdiff(fixed_equations, eq)
    if (length(bad)) stop(fn, ": unknown equation name(s) in 'fixed_equations': ", paste(bad, collapse = ", "), ".", call. = FALSE)
    e <- fixed_equations
  } else {
    i <- as.integer(fixed_equations)
    if (anyNA(i) || any(i < 1) || any(i > length(eq))) stop(fn, ": 'fixed_equations' must give equation indices between 1 and ", length(eq), " or names.", call. = FALSE)
    e <- eq[i]
  }
  which(sub(":.*$", "", coef_names) %in% e)
}

#' @title Restriction matrix for a partial structural change model
#' @description Builds \code{R} in \code{beta = R delta} when only the
#' coefficients in \code{changing} differ across regimes: \code{delta} holds the
#' \code{p} coefficients of the first regime followed by the changing
#' coefficients of each later regime.
#' @param p number of coefficients per regime.
#' @param m number of breaks.
#' @param changing indices of the coefficients allowed to change.
#' @return Matrix with \code{p (m + 1)} rows and \code{p + m * length(changing)} columns.
#' @keywords internal
restriction_matrix <- function(p, m, changing){
  k <- length(changing)
  R <- matrix(0, p * (m + 1), p + k * m)
  R[seq_len(p), seq_len(p)] <- diag(p)
  for (j in seq_len(m) + 1){
    for (cc in seq_len(p)){
      pos <- match(cc, changing)
      if (is.na(pos)) R[(j - 1) * p + cc, cc] <- 1 else R[(j - 1) * p + cc, p + (j - 2) * k + pos] <- 1
    }
  }
  R
}
