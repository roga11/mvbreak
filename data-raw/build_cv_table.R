# Builds the critical-value table shipped in R/sysdata.rda (item B7).
#
# For every number of changing parameters q = 1..q_max and every trimming in
# `trms`, one simulation of B paths with n steps gives the sup-LR draws for
# 1 to M breaks on common paths, from which the supLR, UDmax, WDmax and
# sequential critical values are taken (see critical_values()). The largest
# number of breaks per trimming follows Bai and Perron (2003). Each (q, trm)
# cell is simulated sequentially from its own seed, so the table does not
# depend on how many workers are used.
#
# Usage, from the package directory with the package installed:
#   Rscript data-raw/build_cv_table.R [workers]
# Takes about 10 CPU-hours at the default settings.

suppressMessages(library(mvbreak))
attach(asNamespace("mvbreak"), name = "mvbreak_internals", warn.conflicts = FALSE)   # the ported internals are not exported
args <- commandArgs(trailingOnly = TRUE)
workers <- if (length(args) >= 1) as.integer(args[1]) else 1L

q_max  <- 60L
trms   <- c(0.05, 0.10, 0.15, 0.20, 0.25)
m_max  <- c(9L, 8L, 5L, 3L, 2L)
levels <- c(0.90, 0.95, 0.975, 0.99)
B <- 10000L; n <- 1000L
base_seed <- 250493L

cells <- expand.grid(q = seq_len(q_max), it = seq_along(trms))
one_cell <- function(k){
  q <- cells$q[k]; it <- cells$it[k]; M <- m_max[it]
  h <- mvbreak:::trim_length(trms[it], n)
  set.seed(base_seed + 1000L * it + q)
  draws <- mvbreak:::sim_supF_draws(M, q, h, n, B)
  # supLR and Seq for 1..M breaks come from one call; UDmax and WDmax depend on
  # the largest number of breaks allowed, so they are computed for every maximum
  # from the same draws
  full <- mvbreak:::cv_from_draws(draws, q, c("supLR", "seq"), levels)
  list(supLR = full[, paste0("supLR(m=", 1:M, ")"), drop = FALSE],
       seq   = full[, paste0("Seq(", 2:(M + 1), "|", 1:M, ")"), drop = FALSE],
       UDmax = matrix(unlist(lapply(seq_len(M), function(Mx) cv_from_draws_col(draws, q, Mx, "UDmax", levels))), length(levels), M),
       WDmax = matrix(unlist(lapply(seq_len(M), function(Mx) cv_from_draws_col(draws, q, Mx, "WDmax", levels))), length(levels), M))
}
cv_from_draws_col <- function(draws, q, Mx, test, levels){
  mvbreak:::cv_from_draws(draws[, 1:Mx, drop = FALSE], q, test, levels)[, test]
}
t0 <- proc.time()[3]
if (workers > 1){
  cl <- parallel::makePSOCKcluster(workers)
  parallel::clusterExport(cl, c("cells", "m_max", "trms", "levels", "B", "n", "base_seed", "cv_from_draws_col"))
  res <- parallel::parLapply(cl, seq_len(nrow(cells)), one_cell)
  parallel::stopCluster(cl)
} else {
  res <- lapply(seq_len(nrow(cells)), one_cell)
}
cat(sprintf("simulated %d cells in %.1f minutes\n", nrow(cells), (proc.time()[3] - t0) / 60))

M_all <- max(m_max)
arr <- function() array(NA_real_, dim = c(q_max, length(trms), M_all, length(levels)))
.cv_table <- list(trm = trms, level = levels, q_max = q_max, m_max = m_max, B = B, n = n,
                  supLR = arr(), UDmax = arr(), WDmax = arr(), seq = arr())
for (k in seq_len(nrow(cells))){
  q <- cells$q[k]; it <- cells$it[k]; M <- m_max[it]; r <- res[[k]]
  for (mm in seq_len(M)){
    .cv_table$supLR[q, it, mm, ] <- r$supLR[, mm]
    .cv_table$seq[q, it, mm, ]   <- r$seq[, mm]
    .cv_table$UDmax[q, it, mm, ] <- r$UDmax[, mm]
    .cv_table$WDmax[q, it, mm, ] <- r$WDmax[, mm]
  }
}
save(.cv_table, file = "R/sysdata.rda", compress = "xz")
cat("written R/sysdata.rda\n")
