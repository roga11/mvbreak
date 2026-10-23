# S3 methods for the fitted models (classes "mvbreakEst" and "mvbreakOrdered") and
# summaries of the test objects ("mvbreakTest", "mvbreakSelect").

# ---- helpers ------------------------------------------------------------------

#' @title Regime of each observation
#' @description Index of the regime, 1 to \code{m + 1}, each observation belongs to.
#' @param dates break dates (last observation of each regime but the last).
#' @param bigt number of observations.
#' @return Integer vector of length \code{bigt}.
#' @keywords internal
regime_index <- function(dates, bigt){
  findInterval(seq_len(bigt), c(0, dates), left.open = TRUE)
}

#' @title Coefficient table with standard errors
#' @description Builds the estimate / std. error / t value / p value table of
#' \code{summary()} from a coefficient vector and its covariance matrix.
#' @param est coefficient estimates.
#' @param V covariance matrix of the estimates.
#' @return Matrix with four columns.
#' @keywords internal
coef_table <- function(est, V){
  se <- sqrt(pmax(diag(V), 0))
  tval <- est / se
  tab <- cbind(Estimate = est, `Std. Error` = se, `t value` = tval, `Pr(>|t|)` = 2 * stats::pnorm(-abs(tval)))
  rownames(tab) <- names(est)
  tab
}

#' @title Shade an interval and mark a break date on a plot
#' @keywords internal
#' @param date break date.
#' @param bounds optional lower and upper bound.
#' @param col colour of the shading.
mark_break <- function(date, bounds = NULL, col = grDevices::adjustcolor("steelblue", 0.25)){
  if (!is.null(bounds) && all(is.finite(bounds))){
    usr <- graphics::par("usr")
    graphics::rect(bounds[1], usr[3], bounds[2], usr[4], col = col, border = NA)
  }
  graphics::abline(v = date, lty = 2, col = "steelblue4")
}

# ---- class "mvbreakEst" --------------------------------------------------------

#' @title Methods for fitted structural change models
#'
#' @description Standard methods for the objects returned by
#' \code{breaks_varp()}, \code{breaks_mvreg()}, \code{breaks_sur()},
#' \code{breaks_panel()} and \code{breaks_system()}.
#'
#' @details \code{coef()} returns the \code{p} by \code{m + 1} matrix of regime
#' coefficients. \code{vcov()} returns the asymptotic covariance matrix of
#' the stacked coefficients of all regimes given the estimated break dates,
#' \code{R solve(t(R) t(Zbar) solve(Omega) Zbar R) t(R)}, with the regime covariance
#' matrices in \code{Omega}; because the break dates converge faster than the
#' coefficients, this is the covariance matrix of Theorem 2 of Qu and Perron
#' (2007). \code{summary()} adds standard errors and t statistics for every
#' regime. \code{residuals()} and \code{fitted()} are \code{T} by \code{n}
#' matrices on the effective sample (for a VAR, after the initial values).
#' \code{logLik()} is the Gaussian log-likelihood at the estimates, with the
#' number of free coefficients, the break dates and the covariance parameters
#' as degrees of freedom. \code{confint()} gives intervals for the break dates
#' (\code{parm = "dates"}, through \code{confint_breaks()}) or normal intervals
#' for the coefficients (\code{parm = "coef"}). \code{plot()} draws each
#' dependent variable with the break dates as dashed lines and, if \code{ci} is
#' the output of \code{confint_breaks()}, shades the intervals of its first level.
#'
#' @param x,object an object of class \code{"mvbreakEst"}.
#' @param digits number of significant digits to print.
#' @param parm \code{"dates"} or \code{"coef"}.
#' @param level confidence level(s).
#' @param which which equations to plot, by index or name.
#' @param ci optional matrix from \code{confint_breaks()} to shade.
#' @param ... for \code{confint()} with \code{parm = "dates"}, passed to
#' \code{confint_breaks()}; for \code{plot()}, passed to \code{plot()}; unused otherwise.
#'
#' @return \code{print()} and \code{plot()} return \code{x} invisibly;
#' \code{summary()} an object of class \code{"summary.mvbreakEst"}; the others
#' the usual objects.
#'
#' @references Qu, Z. and P. Perron (2007), "Estimating and Testing Structural Changes in Multivariate Regressions", \emph{Econometrica}, 75, 459-502.
#'
#' @examples
#' y <- as.matrix(data1[, 1:2]); z <- as.matrix(data1[, 3:5])
#' fit <- breaks_sur(y, x = list(z[, 2, drop = FALSE], z[, 3, drop = FALSE]), m = 2, trim = 0.2)
#' fit
#' summary(fit)
#' confint(fit)
#' plot(fit, ci = confint_breaks(fit))
#'
#' @name mvbreakEst-methods
#' @importFrom stats coef vcov residuals fitted nobs logLik confint
NULL

#' @rdname mvbreakEst-methods
#' @export
print.mvbreakEst <- function(x, digits = 4, ...){
  st <- x[["settings"]]; md <- x[["model"]]
  what <- switch(st[["change"]], coef = "coefficients", cov = "covariance matrix", both = "coefficients and covariance matrix")
  cat("Structural changes in a ", md[["type"]], if (!is.null(md[["lags"]])) paste0("(", md[["lags"]], ")"), " with ", md[["n"]], " equations\n", sep = "")
  cat("Change in the ", what, "; ", st[["m"]], " break(s); trimming ", st[["trim"]], " (h = ", st[["h"]], ")",
      if (length(st[["fixed"]])) paste0("; ", length(st[["fixed"]]), " coefficient(s) held fixed") else "", "\n", sep = "")
  cat("\nBreak dates: ", paste(x[["dates"]], collapse = ", "), "\n", sep = "")
  cat("\nCoefficients:\n"); print(round(x[["coef"]], digits))
  cat("\nCovariance matrix of the errors")
  if (st[["change"]] == "coef"){ cat(":\n"); print(round(x[["sigma"]][[1]], digits)) }
  else { cat(", by regime:\n"); for (nm in names(x[["sigma"]])){ cat(nm, "\n"); print(round(x[["sigma"]][[nm]], digits)) } }
  invisible(x)
}

#' @rdname mvbreakEst-methods
#' @export
coef.mvbreakEst <- function(object, ...) object[["coef"]]

#' @rdname mvbreakEst-methods
#' @export
vcov.mvbreakEst <- function(object, ...){
  eng <- object[["engine"]]; n <- eng[["n"]]; bigt <- eng[["bigt"]]; m <- eng[["m"]]; R <- eng[["R"]]
  p <- object[["model"]][["p"]]
  Zbar <- pzbar(eng[["matz"]], m, as.matrix(eng[["dx"]]), bigt) %*% R
  reg <- regime_index(as.integer(eng[["dx"]]), bigt)
  isig <- lapply(seq_len(m + 1), function(j) solve(eng[["rvv"]][((j - 1) * n + 1):(j * n), , drop = FALSE]))
  A <- matrix(0, ncol(R), ncol(R))
  for (t in seq_len(bigt)){
    Zt <- Zbar[((t - 1) * n + 1):(t * n), , drop = FALSE]
    A <- A + t(Zt) %*% isig[[reg[t]]] %*% Zt
  }
  V <- R %*% solve(A) %*% t(R)
  nm <- as.vector(outer(object[["model"]][["coef_names"]], paste0("regime ", seq_len(m + 1)), function(a, b) paste(b, a, sep = ": ")))
  dimnames(V) <- list(nm, nm)
  V
}

#' @rdname mvbreakEst-methods
#' @export
residuals.mvbreakEst <- function(object, ...){
  eng <- object[["engine"]]; n <- eng[["n"]]; bigt <- eng[["bigt"]]
  res <- eng[["maty"]] - pzbar(eng[["matz"]], eng[["m"]], as.matrix(eng[["dx"]]), bigt) %*% eng[["rbeta"]]
  out <- matrix(res, bigt, n, byrow = TRUE); colnames(out) <- object[["model"]][["equations"]]
  out
}

#' @rdname mvbreakEst-methods
#' @export
fitted.mvbreakEst <- function(object, ...){
  y <- object[["data"]][["y"]]; colnames(y) <- object[["model"]][["equations"]]
  y - residuals.mvbreakEst(object)
}

#' @rdname mvbreakEst-methods
#' @export
nobs.mvbreakEst <- function(object, ...) object[["engine"]][["bigt"]]

#' @rdname mvbreakEst-methods
#' @export
logLik.mvbreakEst <- function(object, ...){
  eng <- object[["engine"]]; n <- eng[["n"]]; bigt <- eng[["bigt"]]; m <- eng[["m"]]
  lens <- diff(c(0, as.integer(eng[["dx"]]), bigt))
  ll <- 0
  for (j in seq_len(m + 1)){
    S <- eng[["rvv"]][((j - 1) * n + 1):(j * n), , drop = FALSE]
    ll <- ll - lens[j] * n / 2 * (log(2 * pi) + 1) - lens[j] / 2 * log(det(S))
  }
  ncov <- n * (n + 1) / 2 * if (eng[["brv"]] == 1) m + 1 else 1
  structure(ll, df = ncol(eng[["R"]]) + m + ncov, nobs = bigt, class = "logLik")
}

#' @rdname mvbreakEst-methods
#' @export
confint.mvbreakEst <- function(object, parm = c("dates", "coef"), level = 0.95, ...){
  parm <- match.arg(parm)
  if (parm == "dates") return(confint_breaks(object, level = level, ...))
  est <- as.vector(object[["coef"]]); V <- vcov.mvbreakEst(object); se <- sqrt(pmax(diag(V), 0))
  out <- NULL
  for (l in level){
    z <- stats::qnorm(1 - (1 - l) / 2)
    cols <- cbind(est - z * se, est + z * se); colnames(cols) <- paste0(c("lower ", "upper "), 100 * l, "%")
    out <- cbind(out, cols)
  }
  rownames(out) <- rownames(V)
  out
}

#' @rdname mvbreakEst-methods
#' @export
summary.mvbreakEst <- function(object, ...){
  V <- vcov.mvbreakEst(object)
  est <- as.vector(object[["coef"]]); names(est) <- rownames(V)
  structure(list(fit = object, coefficients = coef_table(est, V), logLik = logLik.mvbreakEst(object)),
            class = "summary.mvbreakEst")
}

#' @rdname mvbreakEst-methods
#' @export
print.summary.mvbreakEst <- function(x, digits = 4, ...){
  fit <- x[["fit"]]; m <- fit[["settings"]][["m"]]; p <- fit[["model"]][["p"]]
  print.mvbreakEst(fit, digits = digits)
  cat("\nCoefficients with standard errors (given the break dates):\n")
  tab <- x[["coefficients"]]
  for (j in seq_len(m + 1)){
    cat("regime ", j, "\n", sep = "")
    rows <- tab[((j - 1) * p + 1):(j * p), , drop = FALSE]; rownames(rows) <- fit[["model"]][["coef_names"]]
    stats::printCoefmat(rows, digits = digits, signif.stars = FALSE)
  }
  ll <- x[["logLik"]]
  cat("\nLog-likelihood: ", format(as.numeric(ll), digits = digits + 2), " on ", attr(ll, "df"), " parameters; T = ", attr(ll, "nobs"), "\n", sep = "")
  invisible(x)
}

#' @rdname mvbreakEst-methods
#' @export
plot.mvbreakEst <- function(x, which = NULL, ci = NULL, ...){
  y <- x[["data"]][["y"]]; eq <- x[["model"]][["equations"]]; offset <- if (is.null(x[["engine"]][["offset"]])) 0L else x[["engine"]][["offset"]]
  idx <- if (is.null(which)) seq_along(eq) else if (is.character(which)) match(which, eq) else as.integer(which)
  if (anyNA(idx) || any(idx < 1) || any(idx > length(eq))) stop("plot: 'which' must name equations of the model.", call. = FALSE)
  old <- graphics::par(mfrow = c(length(idx), 1), mar = c(3, 4, 2, 1)); on.exit(graphics::par(old))
  tt <- seq_len(nrow(y)) + offset
  for (i in idx){
    graphics::plot(tt, y[, i], type = "l", xlab = "", ylab = eq[i], main = eq[i], ...)
    for (j in seq_along(x[["dates"]])){
      bounds <- if (!is.null(ci)) ci[j, 2:3] else NULL
      mark_break(x[["dates"]][j], bounds)
    }
  }
  invisible(x)
}

# ---- class "mvbreakOrdered" ----------------------------------------------------

#' @title Methods for locally ordered breaks models
#'
#' @description Standard methods for the objects returned by
#' \code{breaks_ordered()}: \code{print()}, \code{summary()}, \code{coef()},
#' \code{vcov()}, \code{residuals()}, \code{fitted()}, \code{nobs()},
#' \code{logLik()}, \code{confint()} and \code{plot()}, with the same meaning
#' as for \code{"mvbreakEst"} objects (see \code{\link{mvbreakEst-methods}}).
#' \code{coef()} returns the \code{p} by 2 matrix of each coefficient before
#' and after its own equation's break; \code{vcov()} the covariance matrix of
#' those \code{2p} values given the dates; \code{plot()} marks the two dates.
#'
#' @inheritParams mvbreakEst-methods
#' @param x,object an object of class \code{"mvbreakOrdered"}.
#'
#' @return As for \code{\link{mvbreakEst-methods}}.
#'
#' @examples
#' y <- as.matrix(data2[, 1:2]); z <- as.matrix(data2[, 3:5])
#' fit <- breaks_ordered(y, x = list(z[, 2, drop = FALSE], z[, 3, drop = FALSE]), trim = 0.2)
#' fit
#' summary(fit)
#' plot(fit)
#'
#' @name mvbreakOrdered-methods
NULL

#' @rdname mvbreakOrdered-methods
#' @export
print.mvbreakOrdered <- function(x, digits = 4, ...){
  md <- x[["model"]]; st <- x[["settings"]]
  cat("Locally ordered breaks in a ", md[["design"]], if (md[["lags"]] > 0) paste0("(", md[["lags"]], ")"), " with ", md[["n"]], " equations\n", sep = "")
  cat("Equations breaking first: ", paste(md[["first"]], collapse = ", "), "; trimming ", st[["trim"]], " (h = ", st[["h"]], "); max_gap = ", st[["max_gap"]], "\n", sep = "")
  cat("\nBreak dates: k1 = ", x[["dates"]][["k1"]], " (first group), k2 = ", x[["dates"]][["k2"]], " (second group)\n", sep = "")
  cat("LR test of no break against ordered breaks: ", format(x[["lr"]][["statistic"]], digits = digits), " (see test_breaks())\n", sep = "")
  cat("\nCoefficients before and after the break of their own equation:\n"); print(round(x[["coef"]], digits))
  cat("\nCovariance matrix of the errors:\n"); print(round(x[["sigma"]], digits))
  invisible(x)
}

#' @rdname mvbreakOrdered-methods
#' @export
coef.mvbreakOrdered <- function(object, ...) object[["coef"]]

#' @rdname mvbreakOrdered-methods
#' @export
vcov.mvbreakOrdered <- function(object, ...){
  eng <- object[["engine"]]; n <- eng[["n"]]; bigt <- eng[["bigt"]]; p <- object[["model"]][["p"]]
  matz <- eng[["matz"]]; grp1 <- eng[["grp1"]]; k1 <- eng[["k1"]]; k2 <- eng[["k2"]]
  isig <- solve(object[["sigma"]])
  A <- matrix(0, 2 * p, 2 * p)
  for (t in seq_len(bigt)){
    Zt <- matz[((t - 1) * n + 1):(t * n), , drop = FALSE]
    shift <- ifelse(grp1, t > k1, t > k2)
    Xt <- cbind(Zt, sweep(Zt, 2, as.numeric(shift), "*"))
    A <- A + t(Xt) %*% isig %*% Xt
  }
  J <- rbind(cbind(diag(p), matrix(0, p, p)), cbind(diag(p), diag(p)))   # (base, shift) -> (before, after)
  V <- J %*% solve(A) %*% t(J)
  nm <- c(paste("before:", object[["model"]][["coef_names"]]), paste("after:", object[["model"]][["coef_names"]]))
  dimnames(V) <- list(nm, nm)
  V
}

#' @rdname mvbreakOrdered-methods
#' @export
residuals.mvbreakOrdered <- function(object, ...){
  out <- object[["engine"]][["res"]]; colnames(out) <- object[["model"]][["equations"]]; out
}

#' @rdname mvbreakOrdered-methods
#' @export
fitted.mvbreakOrdered <- function(object, ...) object[["data"]][["y"]] - residuals.mvbreakOrdered(object)

#' @rdname mvbreakOrdered-methods
#' @export
nobs.mvbreakOrdered <- function(object, ...) object[["engine"]][["bigt"]]

#' @rdname mvbreakOrdered-methods
#' @export
logLik.mvbreakOrdered <- function(object, ...){
  eng <- object[["engine"]]; n <- eng[["n"]]; bigt <- eng[["bigt"]]; p <- object[["model"]][["p"]]
  ll <- -bigt * n / 2 * (log(2 * pi) + 1) - bigt / 2 * log(det(object[["sigma"]]))
  structure(ll, df = 2 * p + 2 + n * (n + 1) / 2, nobs = bigt, class = "logLik")
}

#' @rdname mvbreakOrdered-methods
#' @export
confint.mvbreakOrdered <- function(object, parm = c("dates", "coef"), level = 0.95, ...){
  parm <- match.arg(parm)
  if (parm == "dates") return(confint_breaks(object, level = level, ...))
  est <- as.vector(object[["coef"]]); V <- vcov.mvbreakOrdered(object); se <- sqrt(pmax(diag(V), 0))
  out <- NULL
  for (l in level){
    z <- stats::qnorm(1 - (1 - l) / 2)
    cols <- cbind(est - z * se, est + z * se); colnames(cols) <- paste0(c("lower ", "upper "), 100 * l, "%")
    out <- cbind(out, cols)
  }
  rownames(out) <- rownames(V)
  out
}

#' @rdname mvbreakOrdered-methods
#' @export
summary.mvbreakOrdered <- function(object, ...){
  V <- vcov.mvbreakOrdered(object)
  est <- as.vector(object[["coef"]]); names(est) <- rownames(V)
  structure(list(fit = object, coefficients = coef_table(est, V), logLik = logLik.mvbreakOrdered(object)),
            class = "summary.mvbreakOrdered")
}

#' @rdname mvbreakOrdered-methods
#' @export
print.summary.mvbreakOrdered <- function(x, digits = 4, ...){
  fit <- x[["fit"]]; p <- fit[["model"]][["p"]]
  print.mvbreakOrdered(fit, digits = digits)
  cat("\nCoefficients with standard errors (given the break dates):\n")
  tab <- x[["coefficients"]]
  for (j in 1:2){
    cat(c("before", "after")[j], "\n", sep = "")
    rows <- tab[((j - 1) * p + 1):(j * p), , drop = FALSE]; rownames(rows) <- fit[["model"]][["coef_names"]]
    stats::printCoefmat(rows, digits = digits, signif.stars = FALSE)
  }
  ll <- x[["logLik"]]
  cat("\nLog-likelihood: ", format(as.numeric(ll), digits = digits + 2), " on ", attr(ll, "df"), " parameters; T = ", attr(ll, "nobs"), "\n", sep = "")
  invisible(x)
}

#' @rdname mvbreakOrdered-methods
#' @export
plot.mvbreakOrdered <- function(x, which = NULL, ci = NULL, ...){
  y <- x[["data"]][["y"]]; eq <- x[["model"]][["equations"]]; offset <- x[["engine"]][["offset"]]
  idx <- if (is.null(which)) seq_along(eq) else if (is.character(which)) match(which, eq) else as.integer(which)
  if (anyNA(idx) || any(idx < 1) || any(idx > length(eq))) stop("plot: 'which' must name equations of the model.", call. = FALSE)
  old <- graphics::par(mfrow = c(length(idx), 1), mar = c(3, 4, 2, 1)); on.exit(graphics::par(old))
  tt <- seq_len(nrow(y)) + offset
  first_eq <- eq %in% x[["model"]][["first"]]
  for (i in idx){
    graphics::plot(tt, y[, i], type = "l", xlab = "", ylab = eq[i], main = paste0(eq[i], if (first_eq[i]) " (breaks first)" else " (breaks second)"), ...)
    j <- if (first_eq[i]) 1 else 2
    mark_break(x[["dates"]][j], if (!is.null(ci)) ci[j, 2:3] else NULL)
  }
  invisible(x)
}

# ---- tests and selection ------------------------------------------------------

#' @title Summaries of the test and selection objects
#'
#' @description \code{summary()} of a \code{"mvbreakTest"} object returns one
#' table with every statistic, its critical values and the strictest level at
#' which it rejects; for a \code{"mvbreakSelect"} object the same table with
#' the chosen number of breaks as an attribute.
#'
#' @param object an object of class \code{"mvbreakTest"} or \code{"mvbreakSelect"}.
#' @param ... unused.
#'
#' @return A data frame with one row per test.
#'
#' @name mvbreakTest-methods
NULL

#' @rdname mvbreakTest-methods
#' @export
summary.mvbreakTest <- function(object, ...){
  lev <- object[["settings"]][["level"]]
  rows <- list()
  add <- function(name, stat, cvs){
    r <- which(stat > cvs)
    rows[[length(rows) + 1]] <<- data.frame(test = name, statistic = stat, t(matrix(cvs, dimnames = list(paste0("cv ", 100 * lev, "%"), NULL))),
                                            reject_at = if (!is.finite(stat)) NA_character_ else if (length(r)) paste0(100 * (1 - max(lev[r])), "%") else "none",
                                            check.names = FALSE, stringsAsFactors = FALSE)
  }
  if (!is.null(object[["supLR"]])) for (i in seq_len(nrow(object[["supLR"]]))) add(paste0(object[["test"]], " ", rownames(object[["supLR"]])[i]), object[["supLR"]][i, 1], object[["supLR"]][i, -1])
  if (!is.null(object[["UDmax"]])) add("UDmax", object[["UDmax"]][1, 1], object[["UDmax"]][1, -1])
  if (!is.null(object[["WDmax"]])) for (i in seq_len(nrow(object[["WDmax"]]))){
    cvs <- rep(NA_real_, length(lev)); cvs[i] <- object[["WDmax"]][i, "cv"]
    add(paste0("WDmax at ", rownames(object[["WDmax"]])[i]), object[["WDmax"]][i, "statistic"], cvs)
  }
  if (!is.null(object[["seq"]])) for (i in seq_len(nrow(object[["seq"]]))) add(rownames(object[["seq"]])[i], object[["seq"]][i, 1], object[["seq"]][i, -1])
  out <- do.call(rbind, rows); rownames(out) <- NULL
  out
}

#' @rdname mvbreakTest-methods
#' @export
summary.mvbreakSelect <- function(object, ...){
  out <- summary.mvbreakTest(object[["tests"]])
  attr(out, "n_breaks") <- object[["n_breaks"]]
  out
}
