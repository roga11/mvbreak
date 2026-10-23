# Locally ordered breaks: the reference run of obreak2 / brorder under Octave on both
# data sets, and the argmax used for the intervals against a brute-force search.

test_that("breaks_ordered reproduces the reference on data2 (dates 43, 53; LR 58.019)", {
  y <- as.matrix(data2[, 1:2]); z <- as.matrix(data2[, 3:5])
  fit <- breaks_ordered(y, x = list(z[, 2, drop = FALSE], z[, 3, drop = FALSE]), trim = 0.2)
  expect_s3_class(fit, "mvbreakOrdered")
  expect_equal(unname(fit$dates), c(43L, 53L))
  expect_equal(fit$lr$statistic, 58.018892540059923, tolerance = 1e-8)
  expect_equal(fit$settings$max_gap, 10L)
  expect_equal(dim(fit$coef), c(4, 2))
  tb <- test_breaks(fit)
  expect_equal(unname(tb$supLR[1, "statistic"]), fit$lr$statistic)
  expect_equal(tb$q, 4L)
})

test_that("breaks_ordered reproduces the reference on data1 (dates 38, 56; LR 19.532)", {
  y <- as.matrix(data1[, 1:2]); z <- as.matrix(data1[, 3:5])
  fit <- breaks_ordered(y, x = list(z[, 2, drop = FALSE], z[, 3, drop = FALSE]), trim = 0.2)
  expect_equal(unname(fit$dates), c(38L, 56L))
  expect_equal(fit$lr$statistic, 19.531841186155788, tolerance = 1e-8)
})

test_that("the three designs and 'first' by name work", {
  y <- as.matrix(data2[, 1:2]); colnames(y) <- c("a", "b"); z <- as.matrix(data2[, 3:5])
  f_mv <- breaks_ordered(y, x = z[, 2:3], trim = 0.2)
  expect_equal(f_mv$model$design, "multivariate regression")
  expect_equal(nrow(f_mv$coef), 6)
  f_var <- breaks_ordered(y, lags = 1, first = "b", trim = 0.2)
  expect_equal(f_var$model$design, "VAR")
  expect_equal(f_var$model$first, "b")
  expect_equal(nrow(residuals(f_var)), nrow(y) - 1)
  expect_error(breaks_ordered(y, first = c("a", "b")), "at least one must break second")
  expect_error(breaks_ordered(y[, 1, drop = FALSE]), "at least two equations")
})

test_that("the one-pass argmax equals a brute-force search on the same paths", {
  # small grid so that the brute force is cheap
  M <- 2; step <- 0.1; N <- M / step; v <- seq(-M, M, by = step)
  P <- list(P1j = 3, P2j = 2, P12j = 0.4, P1jp1 = 2.5, P2jp1 = 1.5, P12jp1 = -0.3)
  brute <- function(W){
    sbv1 <- c(rev(W[, 1]), 0, W[, 4]); sbv2 <- c(rev(W[, 2]), 0, W[, 3])
    V1 <- matrix(v, length(v), length(v)); V2 <- t(V1)
    Bv1 <- matrix(sbv1, length(v), length(v)); Bv2 <- t(matrix(sbv2, length(v), length(v)))
    Z <- with(P, (V1 < 0) * (V2 <= 0) * (V1 <= V2) * ((Bv1 + V1 / 2) + sqrt(P2j / P1j) * Bv2 + V2 / 2 * (P2j + 2 * P12j) / P1j) +
                 (V1 < 0) * (V2 > 0) * (Bv1 + V1 / 2 + sqrt(P2jp1 / P1j) * Bv2 - V2 / 2 * (P2jp1 / P1j)) +
                 (V1 >= 0) * (V2 > 0) * (V1 <= V2) * (sqrt(P1jp1 / P1j) * Bv1 + sqrt(P2jp1 / P1j) * Bv2 - (P2jp1 / P1j) * V2 / 2 - ((P1jp1 + 2 * P12jp1) / P1j) * V1 / 2) +
                 (V1 > V2) * (-1))
    i <- which(Z == max(Z), arr.ind = TRUE)[1, ]
    c(V1[i[1], i[2]], V2[i[1], i[2]])
  }
  # the C++ routine draws its own paths; compare distributions of the argmax instead of paths
  set.seed(3)
  bf <- t(replicate(300, brute(apply(matrix(rnorm(N * 4, sd = sqrt(step)), N, 4), 2, cumsum))))
  set.seed(4)
  op <- mvbreak:::sim_ordered_argmax(P$P1j, P$P2j, P$P12j, P$P1jp1, P$P2jp1, P$P12jp1, M, step, 3000L)
  expect_equal(mean(bf[, 1]), mean(op[, 1]), tolerance = 0.25)
  expect_equal(mean(bf[, 2]), mean(op[, 2]), tolerance = 0.25)
  expect_true(all(abs(op) <= M + 1e-12))
})

test_that("intervals for ordered breaks are reproducible with a seed and reject robust options", {
  y <- as.matrix(data2[, 1:2]); z <- as.matrix(data2[, 3:5])
  fit <- breaks_ordered(y, x = list(z[, 2, drop = FALSE], z[, 3, drop = FALSE]), trim = 0.2)
  a <- confint_breaks(fit, B = 2000, seed = 1); b <- confint_breaks(fit, B = 2000, seed = 1)
  expect_identical(a, b)
  expect_equal(rownames(a), c("k1", "k2"))
  expect_true(all(a[, "lower 90%"] <= a[, "estimate"] & a[, "estimate"] <= a[, "upper 90%"]))
  expect_error(confint_breaks(fit, robust = TRUE), "serially uncorrelated")
  expect_error(test_breaks(fit, robust = TRUE), "no HAC-robust test")
  expect_error(select_breaks(fit), "must be a model fitted by")
})
