#' @title Critical values for the structural change tests
#'
#' @description Critical values of the sup-LR test against a fixed number of
#' breaks, of the UDmax and WDmax tests against an unknown number of breaks up
#' to a maximum, and of the sequential test Seq(l+1 | l), for a given number of
#' changing parameters and trimming. Values are read from the table shipped with
#' the package when the case is covered, and simulated otherwise (or always, with
#' \code{method = "simulate"}). The simulation approximates the limiting
#' distributions of Theorems 5 and 6 and Corollary 2 of Qu and Perron (2007) by
#' discretised Wiener processes, as in Bai and Perron (1998, 2003); all four
#' tests are computed from the same simulated paths.
#'
#' @param q number of parameters that change at each break: the coefficients
#' allowed to change plus, when the covariance matrix changes, its
#' \code{n(n+1)/2} distinct elements.
#' @param m number of breaks. For \code{"supLR"} and \code{"seq"}, values for
#' 1 to \code{m} breaks are returned; for \code{"UDmax"} and \code{"WDmax"},
#' \code{m} is the largest number of breaks the test allows.
#' @param trm trimming, as a fraction of the sample.
#' @param test which tests to return critical values for; any subset of
#' \code{"supLR"}, \code{"UDmax"}, \code{"WDmax"}, \code{"seq"}.
#' @param level confidence levels.
#' @param method \code{"table"} (the default) reads the shipped table and
#' simulates only the cases it does not cover; \code{"simulate"} always simulates.
#' @param B number of simulated paths when simulating.
#' @param n number of steps used to approximate the Wiener process when simulating.
#' @param seed optional seed for the simulation. The caller's random number
#' state is left unchanged.
#' @param workers number of workers for the simulation; 1 is sequential.
#' Results are reproducible for a given seed and number of workers.
#'
#' @details All values are on the scale of the likelihood-ratio statistics
#' computed by the package: \code{supLR(m) = 2 (lr_m - lr_0)}. The UDmax
#' statistic is \code{max_m supLR(m) / m} and the WDmax statistic is
#' \code{max_m [c(1) / c(m)] supLR(m)}, where \code{c(m)} is the sup-LR critical
#' value at the same level, so that each number of breaks receives equal
#' weight (Bai and Perron 1998). The Seq(l+1 | l) critical value is the
#' \code{level^(1/(l+1))} quantile of the one-break distribution (Qu and Perron
#' 2007, Theorem 6).
#'
#' @return A matrix with one row per level and one column per requested
#' quantity: \code{supLR(m=1)}, ..., \code{UDmax}, \code{WDmax},
#' \code{Seq(2|1)}, .... Attributes \code{q}, \code{trm}, \code{method},
#' \code{B} and \code{n} record how the values were obtained.
#'
#' @references Bai, J. and P. Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes", \emph{Econometrica}, 66, 47-78.
#' @references Bai, J. and P. Perron (2003), "Critical Values for Multiple Structural Change Tests", \emph{Econometrics Journal}, 6, 72-78.
#' @references Qu, Z. and P. Perron (2007), "Estimating and Testing Structural Changes in Multivariate Regressions", \emph{Econometrica}, 75, 459-502.
#'
#' @examples
#' critical_values(q = 4, m = 2, trm = 0.20, method = "simulate", B = 500, n = 300, seed = 1)
#'
#' @export
critical_values <- function(q, m = 1, trm = 0.15,
                            test = c("supLR", "UDmax", "WDmax", "seq"),
                            level = c(0.90, 0.95, 0.975, 0.99),
                            method = c("table", "simulate"),
                            B = 10000, n = 1000, seed = NULL, workers = 1){
  fn <- "critical_values"
  test <- match.arg(test, several.ok = TRUE)
  method <- match.arg(method)
  if (length(q) != 1 || is.na(q) || q < 1 || q != round(q)) stop(fn, ": 'q' must be a positive integer.", call. = FALSE)
  if (length(m) != 1 || is.na(m) || m < 1 || m != round(m)) stop(fn, ": 'm' must be a positive integer.", call. = FALSE)
  if (length(trm) != 1 || is.na(trm) || trm <= 0 || trm >= 0.5) stop(fn, ": 'trm' must lie strictly between 0 and 0.5.", call. = FALSE)
  if (any(is.na(level)) || any(level <= 0) || any(level >= 1)) stop(fn, ": 'level' must lie strictly between 0 and 1.", call. = FALSE)
  if ((m + 1) * trm >= 1) stop(fn, ": with trimming ", trm, " at most ", ceiling(1 / trm) - 2, " break(s) can be considered; ", m, " requested.", call. = FALSE)

  from_table <- NULL
  if (method == "table"){
    from_table <- cv_from_table(q, m, trm, test, level)
    if (is.null(from_table)){
      message(fn, ": the table does not cover q = ", q, ", m = ", m, ", trm = ", trm, " at these levels; simulating.")
    }
  }
  if (!is.null(from_table)){
    return(from_table)
  }
  h <- trim_length(trm, n)
  if (h < 1 || (m + 1) * h >= n) stop(fn, ": 'n' is too small for this trimming and number of breaks.", call. = FALSE)
  draws <- sim_draws(m, q, h, n, B, seed, workers, fn)
  out <- cv_from_draws(draws, q, test, level)
  attr(out, "q") <- q; attr(out, "trm") <- trm; attr(out, "method") <- "simulate"
  attr(out, "B") <- B; attr(out, "n") <- n
  out
}

#' @title Critical values from simulated draws
#'
#' @description Turns a matrix of simulated sup-LR draws (one column per
#' number of breaks) into the critical values documented in
#' \code{critical_values()}.
#'
#' @param draws \code{B} by \code{M} matrix from \code{sim_supF_draws()}.
#' @param q number of changing parameters.
#' @param test tests requested.
#' @param level confidence levels.
#'
#' @return Matrix, levels by quantities.
#'
#' @keywords internal
cv_from_draws <- function(draws, q, test, level){
  M <- ncol(draws)
  supLR <- apply(draws, 2, stats::quantile, probs = level, names = FALSE, type = 1)
  supLR <- matrix(supLR, length(level), M)
  out <- NULL
  if ("supLR" %in% test){
    cols <- supLR; colnames(cols) <- paste0("supLR(m=", 1:M, ")"); out <- cbind(out, cols)
  }
  if ("UDmax" %in% test){
    ud <- apply(sweep(draws, 2, 1:M, "/"), 1, max)
    cols <- matrix(stats::quantile(ud, level, names = FALSE, type = 1), ncol = 1); colnames(cols) <- "UDmax"; out <- cbind(out, cols)
  }
  if ("WDmax" %in% test){
    wd <- sapply(seq_along(level), function(i){
      w <- supLR[i, 1] / supLR[i, ]
      stats::quantile(apply(sweep(draws, 2, w, "*"), 1, max), level[i], names = FALSE, type = 1)
    })
    cols <- matrix(wd, ncol = 1); colnames(cols) <- "WDmax"; out <- cbind(out, cols)
  }
  if ("seq" %in% test){
    cols <- sapply(1:M, function(l) stats::quantile(draws[, 1], level^(1 / (l + 1)), names = FALSE, type = 1))
    cols <- matrix(cols, length(level), M); colnames(cols) <- paste0("Seq(", 2:(M + 1), "|", 1:M, ")"); out <- cbind(out, cols)
  }
  rownames(out) <- paste0(100 * level, "%")
  out
}

#' @title Critical values from the shipped table
#'
#' @description Looks up the requested critical values in the table built by
#' \code{data-raw/build_cv_table.R}. Returns \code{NULL} when any requested
#' value is not in the table.
#'
#' @inheritParams cv_from_draws
#' @param m number of breaks.
#' @param trm trimming.
#'
#' @return Matrix, levels by quantities, or \code{NULL}.
#'
#' @keywords internal
cv_from_table <- function(q, m, trm, test, level){
  tab <- cv_table_internal()
  if (is.null(tab)) return(NULL)
  it <- match(trm, tab$trm); il <- match(level, tab$level)
  if (is.na(it) || any(is.na(il)) || q > tab$q_max || m > tab$m_max[it]) return(NULL)
  out <- NULL
  if ("supLR" %in% test){
    cols <- tab$supLR[q, it, 1:m, il, drop = TRUE]; cols <- matrix(cols, length(level), m, byrow = TRUE)
    colnames(cols) <- paste0("supLR(m=", 1:m, ")"); out <- cbind(out, cols)
  }
  if ("UDmax" %in% test){
    cols <- matrix(tab$UDmax[q, it, m, il], ncol = 1); colnames(cols) <- "UDmax"; out <- cbind(out, cols)
  }
  if ("WDmax" %in% test){
    cols <- matrix(tab$WDmax[q, it, m, il], ncol = 1); colnames(cols) <- "WDmax"; out <- cbind(out, cols)
  }
  if ("seq" %in% test){
    cols <- tab$seq[q, it, 1:m, il, drop = TRUE]; cols <- matrix(cols, length(level), m, byrow = TRUE)
    colnames(cols) <- paste0("Seq(", 2:(m + 1), "|", 1:m, ")"); out <- cbind(out, cols)
  }
  rownames(out) <- paste0(100 * level, "%")
  attr(out, "q") <- q; attr(out, "trm") <- trm; attr(out, "method") <- "table"
  attr(out, "B") <- tab$B; attr(out, "n") <- tab$n
  out
}

#' @title The shipped critical-value table
#'
#' @description Returns the table stored in \code{R/sysdata.rda}, or
#' \code{NULL} when the package was built without it.
#'
#' @return A list with elements \code{trm}, \code{level}, \code{q_max},
#' \code{m_max}, \code{B}, \code{n} and the arrays \code{supLR}, \code{UDmax},
#' \code{WDmax} and \code{seq}, each indexed by q, trimming, number of breaks and level.
#'
#' @keywords internal
cv_table_internal <- function(){
  if (exists(".cv_table", envir = asNamespace("multiMNbreaks"), inherits = FALSE)){
    get(".cv_table", envir = asNamespace("multiMNbreaks"))
  } else {
    NULL
  }
}

#' @title Simulated sup-LR draws, optionally in parallel and seeded
#'
#' @description Calls \code{sim_supF_draws()} once, or once per worker with
#' the draws split across a socket cluster. When a seed is given the caller's
#' random number state is saved and restored.
#'
#' @inheritParams critical_values
#' @param M largest number of breaks.
#' @param h minimum segment length in steps.
#' @param fn name of the calling function, for messages.
#'
#' @return \code{B} by \code{M} matrix.
#'
#' @keywords internal
sim_draws <- function(M, q, h, n, B, seed, workers, fn){
  workers <- n_workers(workers, fn)
  if (!is.null(seed)){
    had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = globalenv())
    on.exit({
      if (had_seed) assign(".Random.seed", old_seed, envir = globalenv())
      else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
    }, add = TRUE)
  }
  if (workers == 1 || B < 2 * workers){
    if (!is.null(seed)) set.seed(seed)
    return(sim_supF_draws(M, q, h, n, B))
  }
  B_per <- diff(round(seq(0, B, length.out = workers + 1)))
  cl <- parallel::makePSOCKcluster(workers)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterSetRNGStream(cl, iseed = seed)
  pieces <- parallel::clusterMap(cl, function(b, M, q, h, n){
    sim <- utils::getFromNamespace("sim_supF_draws", "multiMNbreaks")
    sim(M, q, h, n, b)
  }, B_per, MoreArgs = list(M = M, q = q, h = h, n = n))
  do.call(rbind, pieces)
}
