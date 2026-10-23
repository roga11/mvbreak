# Tests for the number of breaks on a fitted model, and the sequential procedure
# that selects it. Both take the object returned by the breaks_*() functions and
# read the model, the restrictions and the trimming from it.

#' @title Tests for the number of structural breaks
#'
#' @description Runs the tests of Qu and Perron (2007) on a fitted model: the
#' sup-LR test of no break against \code{k} breaks for \code{k = 1, ...,
#' max_breaks}, the double-maximum tests UDmax and WDmax of no break against an
#' unknown number of breaks up to \code{max_breaks}, and the sequential test
#' Seq(l+1 | l) of \code{l} against \code{l + 1} breaks for \code{l = 1, ...,
#' max_breaks - 1}. Each statistic is reported with its critical values from
#' \code{critical_values()}. The tests allow the same changes (coefficients,
#' covariance matrix or both) and the same restrictions as the fitted model.
#'
#' @param fit a model fitted by \code{breaks_varp()}, \code{breaks_mvreg()},
#' \code{breaks_sur()}, \code{breaks_panel()} or \code{breaks_system()}, or a
#' locally ordered breaks model from \code{breaks_ordered()}, for which only the
#' LR test of no break against ordered breaks exists.
#' @param tests which tests to run: \code{"all"} (the default) or any subset of
#' \code{"supLR"}, \code{"UDmax"}, \code{"WDmax"} and \code{"seq"}.
#' @param max_breaks largest number of breaks considered; by default the number
#' of breaks of the fitted model.
#' @param robust \code{TRUE} to allow serially correlated errors. The sup-LR
#' tests are then replaced by Wald-type sup-F tests with a HAC covariance
#' matrix, which exist only for changes in the coefficients: a model fitted with
#' \code{change = "both"} is tested for changes in its coefficients at the
#' estimated dates, with a message, and a model fitted with \code{change = "cov"}
#' cannot be tested this way.
#' @param prewhit \code{TRUE} to prewhiten the scores with a VAR(1) in the HAC
#' estimate; only used when \code{robust = TRUE}.
#' @param level confidence levels of the critical values.
#' @param method,B,n,seed,workers passed to \code{critical_values()}; \code{workers}
#' is also used for the segment likelihoods.
#'
#' @details The number of parameters that change at each break, \code{q},
#' which indexes the critical values, is the number of coefficients allowed to
#' change plus, when the covariance matrix changes, its \code{n(n+1)/2}
#' distinct elements. For a model with a user-supplied restriction matrix the
#' tests are computed for the fitted number of breaks only, and \code{q}
#' follows the formula of the original code, \code{rows(R)/(m+1) -
#' (rows(R) - cols(R))/m}, which must be an integer.
#'
#' With serially correlated errors the critical values of the tests for changes
#' in the covariance matrix are not valid (Qu and Perron 2007, Section 6.3),
#' which is why \code{robust = TRUE} tests the coefficients only.
#'
#' The sup-LR tests for \code{k} breaks use the break dates estimated for
#' \code{k} breaks; the sequential test Seq(l+1 | l) uses the dates estimated
#' for \code{l} breaks and searches each of the \code{l + 1} segments for one
#' more break with the trimming applied to the segment (Section 6.3).
#'
#' @return An object of class \code{"mvbreakTest"}: a list with
#' \item{supLR}{a \code{max_breaks} by \code{1 + length(level)} matrix: the
#' statistic for each number of breaks and its critical values (sup-F when
#' \code{robust = TRUE});}
#' \item{UDmax}{a 1-row matrix: the statistic and its critical values;}
#' \item{WDmax}{a \code{length(level)} by 2 matrix: the statistic, which depends
#' on the level through its weights, and the critical value at that level;}
#' \item{seq}{a \code{max_breaks - 1} by \code{1 + length(level)} matrix for
#' Seq(2|1), ...; a statistic is \code{NA} when no segment of the partition
#' with \code{l} breaks is long enough to hold another break;}
#' \item{dates}{a list of the break dates estimated for 1, ..., \code{max_breaks} breaks;}
#' \item{q}{the number of changing parameters used for the critical values;}
#' \item{test}{\code{"supLR"} or \code{"supF"};}
#' \item{settings}{\code{max_breaks}, \code{trim}, \code{h}, \code{change},
#' \code{fixed}, \code{robust}, \code{prewhit}, \code{level} and how the critical values were obtained;}
#' \item{call}{the call.}
#' Tests not requested, or not available for the model, are \code{NULL}.
#'
#' @references Bai, J. and P. Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes", \emph{Econometrica}, 66, 47-78.
#' @references Qu, Z. and P. Perron (2007), "Estimating and Testing Structural Changes in Multivariate Regressions", \emph{Econometrica}, 75, 459-502.
#'
#' @seealso \code{select_breaks()} for the sequential procedure that chooses the number of breaks.
#'
#' @examples
#' y <- as.matrix(data1[, 1:2]); z <- as.matrix(data1[, 3:5])
#' fit <- breaks_sur(y, x = list(z[, 2, drop = FALSE], z[, 3, drop = FALSE]), m = 2, trim = 0.2)
#' test_breaks(fit)
#'
#' @export
test_breaks <- function(fit, tests = "all", max_breaks = NULL, robust = FALSE, prewhit = FALSE,
                        level = c(0.90, 0.95, 0.975, 0.99), method = c("table", "simulate"),
                        B = 10000, n = 1000, seed = NULL, workers = 1){
  fn <- "test_breaks"
  if (inherits(fit, "mvbreakOrdered")){
    return(test_ordered(fit, tests, max_breaks, robust, prewhit, level, match.arg(method), B, n, seed, workers, match.call()))
  }
  if (!inherits(fit, "mvbreakEst")) stop(fn, ": 'fit' must be a model fitted by one of the breaks_*() functions.", call. = FALSE)
  all_tests <- c("supLR", "UDmax", "WDmax", "seq")
  if (identical(tests, "all")) tests <- all_tests
  if (!is.character(tests) || !all(tests %in% all_tests)) stop(fn, ": 'tests' must be \"all\" or a subset of ", paste(shQuote(all_tests), collapse = ", "), ".", call. = FALSE)
  for (nm in c("robust", "prewhit")){
    v <- get(nm)
    if (!is.logical(v) || length(v) != 1 || is.na(v)) stop(fn, ": '", nm, "' must be TRUE or FALSE.", call. = FALSE)
  }
  if (any(!is.finite(level)) || any(level <= 0) || any(level >= 1)) stop(fn, ": 'level' must lie strictly between 0 and 1.", call. = FALSE)
  method <- match.arg(method)

  st <- fit[["settings"]]; eng <- fit[["engine"]]; p <- fit[["model"]][["p"]]
  n_eq <- eng[["n"]]; bigt <- eng[["bigt"]]; h <- st[["h"]]; trim <- st[["trim"]]
  change <- st[["change"]]; m_fit <- st[["m"]]
  offset <- if (is.null(eng[["offset"]])) 0L else as.integer(eng[["offset"]])
  if (is.null(max_breaks)) max_breaks <- m_fit
  if (length(max_breaks) != 1 || is.na(max_breaks) || max_breaks < 1 || max_breaks != round(max_breaks)) stop(fn, ": 'max_breaks' must be a positive integer.", call. = FALSE)
  max_breaks <- as.integer(max_breaks)

  if (robust && change == "cov") stop(fn, ": with robust = TRUE only changes in the coefficients can be tested; the model was fitted with change = \"cov\".", call. = FALSE)
  if (robust && change == "both") message(fn, ": with robust = TRUE only the changes in the coefficients are tested, at the dates estimated with change = \"both\".")
  brbeta <- eng[["brbeta"]]; brv <- if (robust) 0L else eng[["brv"]]
  test_name <- if (robust) "supF" else "supLR"

  # restrictions: rebuilt for each number of breaks from 'fixed', or user-supplied (fitted m only)
  R_user <- isTRUE(st[["R_user"]])
  changing <- if (change == "cov") integer(0) else setdiff(seq_len(p), st[["fixed"]])
  if (R_user && max_breaks != m_fit) stop(fn, ": the restriction matrix 'R' of the fit applies to m = ", m_fit, " breaks, so 'max_breaks' must be ", m_fit, ".", call. = FALSE)
  if (R_user){
    unavailable <- intersect(tests, c("UDmax", "WDmax", "seq"))
    if (length(unavailable)){
      if (identical(sort(tests), sort(all_tests))) message(fn, ": with a user-supplied 'R' only the ", test_name, " test for ", m_fit, " break(s) is computed.")
      else stop(fn, ": with a user-supplied 'R' the tests ", paste(unavailable, collapse = ", "), " are not available.", call. = FALSE)
      tests <- "supLR"
    }
  }
  ks <- if (R_user) m_fit else seq_len(max_breaks)
  restriction <- function(k){
    if (R_user) return(st[["R"]])
    if (change == "cov") kronecker(matrix(1, k + 1, 1), diag(p)) else restriction_matrix(p, k, changing)
  }
  unrestricted <- !R_user && (change == "cov" || length(changing) == p)

  # number of changing parameters for the critical values
  q_coef <- if (R_user){
    Rm <- st[["R"]]; nrow(Rm) / (m_fit + 1) - (nrow(Rm) - ncol(Rm)) / m_fit
  } else length(changing)
  q <- q_coef * brbeta + (n_eq * (n_eq + 1) / 2) * brv
  if (abs(q - round(q)) > 1e-8) stop(fn, ": the restriction matrix does not correspond to an integer number of changing parameters per break (", q, ").", call. = FALSE)
  q <- as.integer(round(q))

  floor <- regime_floor(eng[["matz"]], n_eq, prewhit = robust && prewhit)
  check_trimming(h, bigt, max_breaks, floor, fn, warn = FALSE)

  # fits for every number of breaks, sharing the segment likelihoods
  maty <- eng[["maty"]]; matz <- eng[["matz"]]
  bigvec <- if (length(ks) > 1 || ks != m_fit) mlebigvec_par(maty, matz, h, n_eq, bigt, workers) else NULL
  fits <- lapply(ks, function(k){
    if (k == m_fit) return(eng)
    est(maty, matz, n_eq, k, bigt, trim, restriction(k), brbeta, eng[["brv"]], h = h, bigvec = bigvec)
  })
  names(fits) <- ks
  dates <- lapply(fits, function(f) as.integer(f[["dx"]]) + offset)

  cv <- critical_values(q, m = max_breaks, trm = trim, test = all_tests, level = level,
                        method = method, B = B, n = n, seed = seed, workers = workers)
  cv_supLR <- cv[, paste0("supLR(m=", seq_len(max_breaks), ")"), drop = FALSE]
  lev_names <- paste0("cv ", 100 * level, "%")

  stat <- vapply(seq_along(ks), function(i){
    k <- ks[i]; f <- fits[[i]]
    if (!robust) return(suplr(maty, matz, k, n_eq, f[["dx"]], f[["rbeta"]], bigt, brbeta, brv))
    if (unrestricted) as.numeric(supft(maty, matz, k, n_eq, f[["dx"]], f[["rbeta"]], bigt, brbeta, brv, as.integer(prewhit)))
    else as.numeric(suprft(maty, matz, k, n_eq, f[["dx"]], f[["rbeta"]], restriction(k), bigt, brbeta, brv, as.integer(prewhit)))
  }, numeric(1))

  out <- list(supLR = NULL, UDmax = NULL, WDmax = NULL, seq = NULL)
  if ("supLR" %in% tests){
    out[["supLR"]] <- cbind(statistic = stat, t(cv_supLR[, ks, drop = FALSE]))
    dimnames(out[["supLR"]]) <- list(paste0("m=", ks), c("statistic", lev_names))
  }
  if ("UDmax" %in% tests){
    out[["UDmax"]] <- matrix(c(max(stat / ks), cv[, "UDmax"]), 1, dimnames = list("UDmax", c("statistic", lev_names)))
  }
  if ("WDmax" %in% tests){
    wd <- vapply(seq_along(level), function(i) max(stat * cv_supLR[i, 1] / cv_supLR[i, ks]), numeric(1))
    out[["WDmax"]] <- cbind(statistic = wd, cv = cv[, "WDmax"])
    dimnames(out[["WDmax"]]) <- list(paste0(100 * level, "%"), c("statistic", "cv"))
  }
  if ("seq" %in% tests && max_breaks > 1){
    R_ls_1 <- if (unrestricted) NULL else lapply(seq_len(max_breaks), function(i) restriction(1))
    sq <- vapply(seq_len(max_breaks - 1), function(l){
      dt <- fits[[l]][["dx"]]
      if (robust) seqf(maty, matz, n_eq, dt, l + 1, bigt, h, brbeta, brv, as.integer(prewhit), trim, R_ls = R_ls_1)
      else seqlr(maty, matz, n_eq, dt, l + 1, bigt, h, brbeta, brv, trim, R_ls = R_ls_1)
    }, numeric(1))
    cv_seq <- cv[, paste0("Seq(", 2:max_breaks, "|", 1:(max_breaks - 1), ")"), drop = FALSE]
    out[["seq"]] <- cbind(statistic = sq, t(cv_seq))
    dimnames(out[["seq"]]) <- list(paste0("Seq(", 2:max_breaks, "|", 1:(max_breaks - 1), ")"), c("statistic", lev_names))
  }
  out[["dates"]] <- dates
  out[["q"]] <- q
  out[["test"]] <- test_name
  out[["settings"]] <- list(max_breaks = max_breaks, trim = trim, h = h, change = change, fixed = st[["fixed"]],
                            robust = robust, prewhit = prewhit, level = level,
                            critical_values = list(method = attr(cv, "method"), B = attr(cv, "B"), n = attr(cv, "n")))
  out[["call"]] <- match.call()
  class(out) <- "mvbreakTest"
  out
}

#' @title LR test for locally ordered breaks
#' @description The branch of \code{test_breaks()} for a fit of
#' \code{breaks_ordered()}: the LR statistic of no break against ordered
#' breaks with \code{k2 - k1 <= max_gap}, against the critical values of the
#' sup-LR test for one break with all \code{p} coefficients changing.
#' @inheritParams test_breaks
#' @param call the call of \code{test_breaks()}.
#' @return An object of class \code{"mvbreakTest"} with \code{supLR} holding one row.
#' @keywords internal
test_ordered <- function(fit, tests, max_breaks, robust, prewhit, level, method, B, n, seed, workers, call){
  fn <- "test_breaks"
  if (!identical(tests, "all") && !identical(tests, "supLR")) stop(fn, ": only the LR test (\"supLR\") exists for locally ordered breaks.", call. = FALSE)
  if (!is.null(max_breaks) && max_breaks != 1) stop(fn, ": the locally ordered breaks model has one break per group; 'max_breaks' must be 1.", call. = FALSE)
  if (isTRUE(robust)) stop(fn, ": no HAC-robust test exists for locally ordered breaks; the model assumes martingale difference errors.", call. = FALSE)
  if (any(!is.finite(level)) || any(level <= 0) || any(level >= 1)) stop(fn, ": 'level' must lie strictly between 0 and 1.", call. = FALSE)
  st <- fit[["settings"]]; p <- fit[["model"]][["p"]]
  cv <- critical_values(p, m = 1, trm = st[["trim"]], test = "supLR", level = level,
                        method = method, B = B, n = n, seed = seed, workers = workers)
  tab <- matrix(c(fit[["lr"]][["statistic"]], cv[, 1]), 1, dimnames = list("ordered", c("statistic", paste0("cv ", 100 * level, "%"))))
  out <- list(supLR = tab, UDmax = NULL, WDmax = NULL, seq = NULL,
              dates = list(fit[["lr"]][["dates"]]), q = as.integer(p), test = "supLR",
              settings = list(max_breaks = 1L, trim = st[["trim"]], h = st[["h"]], change = "coef", fixed = integer(0),
                              robust = FALSE, prewhit = FALSE, level = level, ordered = TRUE, max_gap = st[["max_gap"]],
                              critical_values = list(method = attr(cv, "method"), B = attr(cv, "B"), n = attr(cv, "n"))),
              call = call)
  class(out) <- "mvbreakTest"
  out
}

#' @title Number of breaks by the sequential procedure
#'
#' @description Chooses the number of breaks as Qu and Perron (2007, Section
#' 6.3) and Bai and Perron (1998) suggest: a double-maximum test (or the sup-LR
#' test for one break) decides whether there is at least one break, and if so
#' the sequential test Seq(l+1 | l) is applied for \code{l = 1, 2, ...} until it
#' no longer rejects or \code{max_breaks} is reached. The model is then
#' refitted with the chosen number of breaks.
#'
#' @inheritParams test_breaks
#' @param level confidence level of the tests, e.g. 0.95 for tests at the 5\% level.
#' @param first_test the test that decides whether there is at least one break:
#' \code{"WDmax"} (the default), \code{"UDmax"} or \code{"supLR"} (the test
#' against one break; the sup-F test when \code{robust = TRUE}).
#'
#' @return An object of class \code{"mvbreakSelect"}: a list with
#' \item{n_breaks}{the chosen number of breaks;}
#' \item{fit}{the model refitted with that number of breaks, or \code{NULL} when it is zero;}
#' \item{tests}{the \code{"mvbreakTest"} object the choice is based on;}
#' \item{level, first_test}{the settings;}
#' \item{call}{the call.}
#' The refit evaluates the call of \code{fit} with \code{m} replaced, so the data
#' must still be available under the names used in that call.
#'
#' @references Bai, J. and P. Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes", \emph{Econometrica}, 66, 47-78.
#' @references Qu, Z. and P. Perron (2007), "Estimating and Testing Structural Changes in Multivariate Regressions", \emph{Econometrica}, 75, 459-502.
#'
#' @examples
#' y <- as.matrix(data1[, 1:2]); z <- as.matrix(data1[, 3:5])
#' fit <- breaks_sur(y, x = list(z[, 2, drop = FALSE], z[, 3, drop = FALSE]), m = 3, trim = 0.2)
#' sel <- select_breaks(fit)
#' sel$n_breaks
#' sel$fit$dates
#'
#' @export
select_breaks <- function(fit, max_breaks = NULL, level = 0.95, first_test = c("WDmax", "UDmax", "supLR"),
                          robust = FALSE, prewhit = FALSE, method = c("table", "simulate"),
                          B = 10000, n = 1000, seed = NULL, workers = 1){
  fn <- "select_breaks"
  if (!inherits(fit, "mvbreakEst")) stop(fn, ": 'fit' must be a model fitted by breaks_varp(), breaks_mvreg(), breaks_sur(), breaks_panel() or breaks_system().", call. = FALSE)
  first_test <- match.arg(first_test)
  if (length(level) != 1 || !is.finite(level) || level <= 0 || level >= 1) stop(fn, ": 'level' must be one number strictly between 0 and 1.", call. = FALSE)
  if (isTRUE(fit[["settings"]][["R_user"]])) stop(fn, ": the sequential procedure needs restrictions given through 'fixed', not a user-supplied 'R'.", call. = FALSE)
  tb <- test_breaks(fit, tests = "all", max_breaks = max_breaks, robust = robust, prewhit = prewhit, level = level,
                    method = match.arg(method), B = B, n = n, seed = seed, workers = workers)
  tb[["call"]] <- NULL
  max_breaks <- tb[["settings"]][["max_breaks"]]
  rejects <- function(row) is.finite(row[["statistic"]]) && row[["statistic"]] > row[[2]]
  at_least_one <- switch(first_test,
                         WDmax = rejects(tb[["WDmax"]][1, ]),
                         UDmax = rejects(tb[["UDmax"]][1, ]),
                         supLR = rejects(tb[["supLR"]][1, ]))
  k <- 0L
  if (at_least_one){
    k <- 1L
    while (k < max_breaks && rejects(tb[["seq"]][k, ])) k <- k + 1L
  }
  refit <- NULL
  if (k == fit[["settings"]][["m"]]){
    refit <- fit
  } else if (k > 0){
    cl <- fit[["call"]]; cl[["m"]] <- k
    refit <- eval(cl, parent.frame())
  }
  structure(list(n_breaks = k, fit = refit, tests = tb, level = level, first_test = first_test, call = match.call()),
            class = "mvbreakSelect")
}

#' @title Print the tests for the number of breaks
#' @description Prints each statistic with its critical values and the
#' strictest level at which it rejects.
#' @param x an object of class \code{"mvbreakTest"}.
#' @param digits number of significant digits.
#' @param ... unused.
#' @return \code{x}, invisibly.
#' @export
print.mvbreakTest <- function(x, digits = 3, ...){
  st <- x[["settings"]]; lev <- st[["level"]]
  what <- switch(st[["change"]], coef = "coefficients", cov = "covariance matrix", both = "coefficients and covariance matrix")
  cat("Tests for structural breaks in the ", what, "\n", sep = "")
  if (st[["robust"]] && st[["change"]] == "both") cat("(HAC-robust tests cover the coefficients only)\n")
  cat("q = ", x[["q"]], " changing parameters per break; trimming ", st[["trim"]], " (h = ", st[["h"]], ")",
      if (st[["robust"]]) paste0("; HAC-robust", if (st[["prewhit"]]) ", prewhitened") else "", "\n", sep = "")
  cv <- st[["critical_values"]]
  cat("Critical values: ", if (identical(cv[["method"]], "table")) "shipped table" else "simulated",
      " (B = ", cv[["B"]], ", n = ", cv[["n"]], ")\n", sep = "")
  reject_at <- function(stat, cvs){
    if (!is.finite(stat)) return("NA")
    r <- which(stat > cvs)
    if (length(r) == 0) "none" else paste0(100 * (1 - max(lev[r])), "%")
  }
  show <- function(title, tab, cv_cols){
    cat("\n", title, "\n", sep = "")
    rej <- vapply(seq_len(nrow(tab)), function(i) reject_at(tab[i, "statistic"], tab[i, cv_cols]), character(1))
    out <- cbind(format(round(tab, digits)), `reject at` = rej)
    print(out, quote = FALSE, right = TRUE)
  }
  if (isTRUE(st[["ordered"]])){
    show(paste0("LR test of no break against locally ordered breaks (k2 - k1 <= ", st[["max_gap"]], ")"), x[["supLR"]], -1)
    return(invisible(x))
  }
  if (!is.null(x[["supLR"]])) show(paste0(x[["test"]], " test of no break against m breaks"), x[["supLR"]], -1)
  if (!is.null(x[["UDmax"]])) show(paste0("UDmax test of no break against up to ", st[["max_breaks"]], " breaks"), x[["UDmax"]], -1)
  if (!is.null(x[["WDmax"]])){
    cat("\nWDmax test of no break against up to ", st[["max_breaks"]], " breaks (weights depend on the level)\n", sep = "")
    tab <- x[["WDmax"]]
    out <- cbind(format(round(tab, digits)), reject = ifelse(tab[, "statistic"] > tab[, "cv"], "yes", "no"))
    print(out, quote = FALSE, right = TRUE)
  }
  if (!is.null(x[["seq"]])) show("Sequential test of l against l+1 breaks", x[["seq"]], -1)
  invisible(x)
}

#' @title Print the selected number of breaks
#' @param x an object of class \code{"mvbreakSelect"}.
#' @param ... passed to \code{print.mvbreakTest()}.
#' @return \code{x}, invisibly.
#' @export
print.mvbreakSelect <- function(x, ...){
  cat("Number of breaks selected: ", x[["n_breaks"]], " (sequential procedure at the ",
      100 * (1 - x[["level"]]), "% level, starting from ", x[["first_test"]], ")\n", sep = "")
  if (!is.null(x[["fit"]])) cat("Break dates: ", paste(x[["fit"]][["dates"]], collapse = ", "), "\n", sep = "")
  cat("\n")
  print(x[["tests"]], ...)
  invisible(x)
}
