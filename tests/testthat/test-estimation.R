# Estimation anchors: the GAUSS output file of mbreak2.prg (Qu and Perron's own run,
# May 2008) and the Matlab translation run under Octave over the validation grid.

y1 <- as.matrix(data1[, 1:2]); z1 <- as.matrix(data1[, 3:5])
sur_x <- list(z1[, 2, drop = FALSE], z1[, 3, drop = FALSE])

test_that("breaks_sur reproduces the GAUSS example (dates 56, 100)", {
  fit <- breaks_sur(y1, sur_x, m = 2, trim = 0.2)
  expect_s3_class(fit, "mvbreakEst")
  expect_equal(fit$dates, c(56L, 100L))
  expect_equal(dim(fit$coef), c(4, 3))
  expect_equal(rownames(fit$coef), c("y1:const", "y1:x1", "y2:const", "y2:x1"))
  expect_equal(length(fit$sigma), 3)
  expect_equal(fit$model$type, "SUR")
})

test_that("breaks_system with an explicit S gives the same fit as breaks_sur", {
  S <- rbind(c(1,0,0,0), c(0,1,0,0), c(0,0,0,0), c(0,0,1,0), c(0,0,0,0), c(0,0,0,1))
  a <- breaks_sur(y1, sur_x, m = 2, trim = 0.2)
  b <- breaks_system(y1, z1, S, m = 2, trim = 0.2)
  expect_equal(b$dates, a$dates)
  expect_equal(unname(b$coef), unname(a$coef))
  expect_equal(unname(b$sigma[[1]]), unname(a$sigma[[1]]))
})

test_that("the restricted model of the GAUSS example gives dates 50, 99", {
  fit <- breaks_sur(y1, sur_x, m = 2, trim = 0.2, fixed = c("y1:x1", "y2:x1"))
  expect_equal(fit$dates, c(50L, 99L))
  # the slopes are the same in every regime
  expect_equal(fit$coef["y1:x1", 1], fit$coef["y1:x1", 2])
  expect_equal(fit$coef["y2:x1", 2], fit$coef["y2:x1", 3])
  expect_equal(fit$settings$fixed, c(2L, 4L))
})

test_that("reference values: SUR with trimming 0.10, two breaks", {
  fit <- breaks_sur(y1, sur_x, m = 2, trim = 0.10)
  expect_equal(fit$dates, c(50L, 100L))
  fit_r <- breaks_sur(y1, sur_x, m = 2, trim = 0.10, fixed = c(2, 4))
  expect_equal(fit_r$dates, c(50L, 99L))
})

test_that("reference values: bivariate VAR(1) designs", {
  # data2, coefficients and covariance change, two breaks, trimming 0.10: dates 52, 62 on the
  # reference's scale, which starts after the initial observation
  fit <- suppressWarnings(breaks_varp(as.matrix(data2[, 1:2]), p = 1, m = 2, trim = 0.10, change = "both"))   # h = 10 < 5 per regressor
  expect_equal(fit$dates, c(52L, 62L) + 1L)
  # data1, covariance only, one break, trimming 0.20: date 55
  fit_cov <- breaks_varp(y1, p = 1, m = 1, trim = 0.20, change = "cov")
  expect_equal(fit_cov$dates, 55L + 1L)
  expect_equal(fit_cov$coef[, 1], fit_cov$coef[, 2])        # coefficients common to both regimes
  expect_false(isTRUE(all.equal(fit_cov$sigma[[1]], fit_cov$sigma[[2]])))
})

test_that("breaks_varp equals breaks_system on hand-built lags", {
  bigt <- nrow(y1)
  z <- cbind(1, y1[1:(bigt - 1), ]); yy <- y1[2:bigt, ]
  a <- breaks_varp(y1, p = 1, m = 2, trim = 0.2)
  b <- breaks_system(yy, z, m = 2, trim = 0.2)
  expect_equal(a$dates, b$dates + 1L)
  expect_equal(unname(a$coef), unname(b$coef))
})

test_that("breaks_mvreg and breaks_panel run and return the documented shapes", {
  f <- breaks_mvreg(y1, z1[, 2:3], m = 1, trim = 0.2)
  expect_equal(dim(f$coef), c(6, 2))
  f0 <- breaks_mvreg(y1, NULL, m = 1, trim = 0.2)                 # means and covariance only
  expect_equal(dim(f0$coef), c(2, 2))
  fp <- breaks_panel(y1, z1[, 2, drop = FALSE], m = 1, trim = 0.2)
  expect_equal(rownames(fp$coef), c("unit1:const", "unit2:const", "x1"))
  fp2 <- breaks_panel(y1, z1[, 2, drop = FALSE], m = 1, trim = 0.2, pooled = FALSE)
  expect_equal(nrow(fp2$coef), 4)
})

test_that("input errors name the function and the argument", {
  expect_error(breaks_sur(y1, sur_x, m = 9, trim = 0.15), "breaks_sur: with T = 150")
  expect_error(breaks_sur(y1, sur_x, m = 2, trim = 0.2, fixed = "y1:nothere"), "unknown coefficient name")
  expect_error(breaks_varp(y1, p = 0), "'p' must be a positive integer")
  expect_error(breaks_mvreg(y1, z1[1:10, ]), "'x' must have as many rows")
  expect_error(breaks_system(y1, z1, m = 2, trim = 0.6), "'trim' must lie")
  expect_message(breaks_sur(y1, sur_x, m = 1, trim = 0.2, change = "cov", fixed = 2), "not used when only the covariance")
})

test_that("h can be passed to reproduce a different rounding of the trimming", {
  # 0.15 * 150 = 22.5: R rounds to 22, GAUSS to 23
  a <- breaks_sur(y1, sur_x, m = 1, trim = 0.15)
  b <- breaks_sur(y1, sur_x, m = 1, trim = 0.15, h = 23)
  expect_equal(a$settings$h, 22L)
  expect_equal(b$settings$h, 23L)
})
