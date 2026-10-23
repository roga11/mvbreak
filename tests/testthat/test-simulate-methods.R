# The simulation functions and the S3 methods.

A1 <- rbind(c(0.5, 0.1), c(0.0, 0.3)); A2 <- rbind(c(-0.2, 0.0), c(0.2, 0.6))

test_that("simulate_varp is reproducible, has the right layout and is recovered by breaks_varp", {
  s <- simulate_varp(200, coef = list(cbind(0, A1), cbind(1, A2)), sigma = diag(2), dates = 100, seed = 1)
  expect_s3_class(s, "mvbreakSim")
  expect_equal(dim(s$y), c(200, 2))
  expect_identical(s$y, simulate_varp(200, coef = list(cbind(0, A1), cbind(1, A2)), sigma = diag(2), dates = 100, seed = 1)$y)
  expect_equal(s$dgp$change, "coef"); expect_equal(s$dgp$m, 1L)
  fit <- breaks_varp(s$y, p = 1, m = 1)
  expect_equal(rownames(s$dgp$coef_matrix), rownames(fit$coef))
  expect_true(abs(fit$dates - 100) <= 2)
  s0 <- simulate_varp(100, coef = cbind(0, A1), sigma = diag(2), seed = 2)
  expect_equal(s0$dgp$m, 0L); expect_null(s0$dgp$dates)
  expect_output(print(s0), "No break")
})

test_that("the other generators match their model functions", {
  s2 <- simulate_mvreg(150, coef = list(rbind(c(0, 1), c(0, -1)), rbind(c(1, 1), c(1, -1))), sigma = list(diag(2), 2 * diag(2)), dates = 75, n_x = 1, seed = 2)
  expect_equal(s2$dgp$change, "both")
  expect_equal(rownames(s2$dgp$coef_matrix), rownames(breaks_mvreg(s2$y, s2$x, m = 1, change = "both")$coef))
  s3 <- simulate_sur(150, coef = list(list(c(0, 1), c(0, 2)), list(c(1, 1), c(1, 2))), sigma = diag(2), dates = 75, n_x = c(1, 1), seed = 3)
  expect_equal(rownames(s3$dgp$coef_matrix), rownames(breaks_sur(s3$y, s3$x, m = 1)$coef))
  s4 <- simulate_panel(150, intercepts = list(c(0, 0, 0), c(1, 1, 1)), slopes = list(0.5, -0.5), sigma = diag(3), dates = 75, n_x = 1, seed = 4)
  expect_equal(rownames(s4$dgp$coef_matrix), rownames(breaks_panel(s4$y, s4$x, m = 1)$coef))
  s5 <- simulate_ordered(120, coef_before = list(c(0, 1), c(0, 1)), coef_after = list(c(1, 1), c(1, 1)), sigma = diag(2), k1 = 50, k2 = 56, n_x = c(1, 1), seed = 5)
  expect_equal(unname(s5$dgp$dates), c(50, 56))
  expect_equal(rownames(s5$dgp$coef_matrix), rownames(breaks_ordered(s5$y, s5$x)$coef))
  expect_error(simulate_varp(100, coef = list(cbind(0, A1)), sigma = diag(2), dates = 50), "m \\+ 1 = 2 elements")
  expect_error(simulate_ordered(100, list(c(0, 1), c(0, 1)), list(c(1, 1), c(1, 1)), diag(2), k1 = 60, k2 = 50), "k1 <= k2")
})

test_that("a seeded generator leaves the caller's RNG state unchanged", {
  set.seed(10); before <- .Random.seed
  invisible(simulate_varp(50, coef = cbind(0, A1), sigma = diag(2), seed = 99))
  expect_identical(.Random.seed, before)
})

y1 <- as.matrix(data1[, 1:2]); z1 <- as.matrix(data1[, 3:5])
sur_x <- list(z1[, 2, drop = FALSE], z1[, 3, drop = FALSE])
fit <- breaks_sur(y1, sur_x, m = 2, trim = 0.2)

test_that("methods on a fitted model are consistent", {
  expect_equal(coef(fit), fit$coef)
  expect_equal(fitted(fit) + residuals(fit), unname(y1), ignore_attr = TRUE)
  expect_equal(nobs(fit), 150L)
  ll <- logLik(fit)
  expect_s3_class(ll, "logLik"); expect_equal(attr(ll, "df"), 17)
  V <- vcov(fit)
  expect_equal(dim(V), c(12, 12)); expect_true(isSymmetric(V, tol = 1e-10))
  s <- summary(fit)
  expect_s3_class(s, "summary.mvbreakEst")
  expect_equal(dim(s$coefficients), c(12, 4))
  expect_output(print(s), "Std. Error")
  expect_output(print(fit), "Break dates: 56, 100")
})

test_that("vcov given the dates equals the GLS covariance of lm on one equation", {
  y <- y1[, 1, drop = FALSE]; f <- breaks_mvreg(y, z1[, 2, drop = FALSE], m = 1, trim = 0.2)
  d <- f$dates; bigt <- nrow(y)
  l1 <- lm(y[1:d] ~ z1[1:d, 2]); l2 <- lm(y[(d + 1):bigt] ~ z1[(d + 1):bigt, 2])
  expect_equal(unname(coef(f)[, 1]), unname(coef(l1)))
  expect_equal(unname(coef(f)[, 2]), unname(coef(l2)))
  s2 <- sum(c(residuals(l1), residuals(l2))^2) / bigt          # the ML estimate of the common variance
  expect_equal(unname(diag(vcov(f))[1:2] / diag(vcov(l1))), rep(s2 / summary(l1)$sigma^2, 2))
})

test_that("plot methods run without error", {
  pdf(NULL); on.exit(dev.off())
  expect_invisible(plot(fit, ci = confint_breaks(fit)))
  expect_invisible(plot(fit, which = "y1"))
  expect_error(plot(fit, which = "nothere"), "must name equations")
  y2 <- as.matrix(data2[, 1:2]); z2 <- as.matrix(data2[, 3:5])
  fo <- breaks_ordered(y2, x = list(z2[, 2, drop = FALSE], z2[, 3, drop = FALSE]), trim = 0.2)
  expect_invisible(plot(fo))
  expect_output(print(summary(fo)), "before")
  expect_equal(attr(logLik(fo), "df"), 13)
})
