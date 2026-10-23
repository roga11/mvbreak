# Confidence intervals for the break dates: the GAUSS output file gives the intervals
# of the unrestricted and restricted models of mbreak2.prg.

y1 <- as.matrix(data1[, 1:2]); z1 <- as.matrix(data1[, 3:5])
sur_x <- list(z1[, 2, drop = FALSE], z1[, 3, drop = FALSE])

test_that("intervals of the GAUSS example are reproduced", {
  fit <- breaks_sur(y1, sur_x, m = 2, trim = 0.2)
  ci <- confint_breaks(fit, level = c(0.95, 0.90))
  expect_equal(unname(ci[, "estimate"]), c(56, 100))
  expect_equal(unname(ci[1, c("lower 95%", "upper 95%")]), c(50, 62))
  expect_equal(unname(ci[2, c("lower 95%", "upper 95%")]), c(94, 106))
  expect_equal(unname(ci[1, c("lower 90%", "upper 90%")]), c(52, 60))
  expect_equal(unname(ci[2, c("lower 90%", "upper 90%")]), c(95, 105))
})

test_that("intervals of the restricted GAUSS example are reproduced", {
  fit <- breaks_sur(y1, sur_x, m = 2, trim = 0.2, fixed = c("y1:x1", "y2:x1"))
  ci <- confint_breaks(fit, level = c(0.95, 0.90))
  expect_equal(unname(ci[1, c("lower 95%", "upper 95%")]), c(43, 57))
  expect_equal(unname(ci[2, c("lower 95%", "upper 95%")]), c(92, 106))
  expect_equal(unname(ci[1, c("lower 90%", "upper 90%")]), c(45, 55))
  expect_equal(unname(ci[2, c("lower 90%", "upper 90%")]), c(94, 104))
})

test_that("robust and prewhitened intervals are finite and contain the estimate", {
  fit <- breaks_sur(y1, sur_x, m = 2, trim = 0.2)
  for (hq in c(FALSE, TRUE)){
    ci <- confint_breaks(fit, robust = TRUE, prewhit = TRUE, hetq = hq)
    expect_true(all(is.finite(ci)))
    expect_true(all(ci[, "lower 90%"] <= ci[, "estimate"] & ci[, "estimate"] <= ci[, "upper 90%"]))
    expect_true(all(ci[, "lower 90%"] >= 1 & ci[, "upper 90%"] <= 150))
  }
})

test_that("VAR intervals are on the original time scale", {
  fit <- breaks_varp(y1, p = 1, m = 1, trim = 0.2)
  ci <- confint_breaks(fit)
  expect_equal(unname(ci[1, "estimate"]), fit$dates)
})

test_that("confint() dispatches to dates and to coefficients", {
  fit <- breaks_sur(y1, sur_x, m = 2, trim = 0.2)
  expect_equal(unname(confint(fit)[, "estimate"]), c(56, 100))
  cc <- confint(fit, "coef", level = 0.9)
  expect_equal(dim(cc), c(12, 2))
  expect_true(all(cc[, 1] < cc[, 2]))
})

test_that("arguments are validated", {
  fit <- breaks_sur(y1, sur_x, m = 2, trim = 0.2)
  expect_error(confint_breaks(fit, level = 1.5), "'level' must lie")
  expect_error(confint_breaks(fit, robust = "yes"), "must be TRUE or FALSE")
  expect_error(confint_breaks(fit, B = 500), "only used for locally ordered breaks")
  expect_warning(confint_breaks(fit, prewhit = TRUE), "no effect")
})

test_that("the quantile search of the limiting distribution is exact and bracket-free", {
  q <- mvbreak:::cvg(1, 1, 1, c(0.025, 0.975))
  expect_equal(q[1], -q[2], tolerance = 1e-6)                  # symmetric case
  # a lopsided case that the fixed bracket of the original code could not solve
  q2 <- mvbreak:::cvg(2, 0.1, 30, c(0.025, 0.975))
  expect_true(all(is.finite(q2)))
  expect_true(is.na(mvbreak:::cvg(-1, 1, 1, 0.5)))
})
