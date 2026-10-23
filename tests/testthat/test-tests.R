# Tests for the number of breaks: anchors from the GAUSS output file and from the
# reference run over the validation grid, and the behaviour of the selection rule.

y1 <- as.matrix(data1[, 1:2]); z1 <- as.matrix(data1[, 3:5])
sur_x <- list(z1[, 2, drop = FALSE], z1[, 3, drop = FALSE])
fit2 <- breaks_sur(y1, sur_x, m = 2, trim = 0.2)

test_that("sup-LR and sequential statistics match the GAUSS output file", {
  tb <- test_breaks(fit2)
  expect_s3_class(tb, "mvbreakTest")
  expect_equal(unname(tb$supLR[, "statistic"]), c(17.939, 55.718), tolerance = 1e-4)
  expect_equal(unname(tb$seq[, "statistic"]), 33.429, tolerance = 1e-4)
  expect_equal(tb$q, 4L)
  expect_equal(tb$test, "supLR")
  # WDmax at 5% with the weights of the response surfaces of the GAUSS code: 34.782
  cv_surface <- sapply(1:2, function(k) mvbreak:::cvlr(k, 4, 0.2)[2])
  wd <- max(tb$supLR[, "statistic"] * cv_surface[1] / cv_surface)
  expect_equal(wd, 34.782, tolerance = 1e-4)
  # UDmax is the largest supLR(m)/m
  expect_equal(unname(tb$UDmax[1, "statistic"]), max(tb$supLR[, "statistic"] / 1:2))
})

test_that("the restricted sup-LR of the GAUSS example is 48.382 with q = 2", {
  fr <- breaks_sur(y1, sur_x, m = 2, trim = 0.2, fixed = c("y1:x1", "y2:x1"))
  tb <- test_breaks(fr)
  expect_equal(unname(tb$supLR[2, "statistic"]), 48.382, tolerance = 1e-4)
  expect_equal(tb$q, 2L)
})

test_that("reference values at trimming 0.10 (sequential test trims within the segment)", {
  # the unfixed code trimmed the sub-segments with the full-sample h and gave other values
  tb <- test_breaks(breaks_sur(y1, sur_x, m = 2, trim = 0.10))
  expect_equal(unname(tb$supLR[, "statistic"]), c(17.278181549050714, 58.034235947555317), tolerance = 1e-8)
  expect_equal(unname(tb$seq[1, "statistic"]), 35.347074967424021, tolerance = 1e-8)
  tv <- test_breaks(suppressWarnings(breaks_varp(as.matrix(data2[, 1:2]), p = 1, m = 2, trim = 0.10, change = "both")))
  expect_equal(unname(tv$supLR[, "statistic"]), c(39.578271295864738, 56.66341867853555), tolerance = 1e-8)
  expect_equal(unname(tv$seq[1, "statistic"]), 33.398739948841921, tolerance = 1e-8)
  expect_equal(tv$q, 9L)                                           # 6 coefficients + 3 covariance elements
})

test_that("the HAC-robust sup-F reproduces the reference", {
  tb <- test_breaks(breaks_sur(y1, sur_x, m = 1, trim = 0.2), robust = TRUE, prewhit = TRUE)
  expect_equal(unname(tb$supLR[1, "statistic"]), 15.333043529289446, tolerance = 1e-8)
  expect_equal(tb$test, "supF")
  expect_error(test_breaks(breaks_sur(y1, sur_x, m = 1, trim = 0.2, change = "cov"), robust = TRUE), "only changes in the coefficients")
  expect_message(test_breaks(breaks_sur(y1, sur_x, m = 1, trim = 0.2, change = "both"), robust = TRUE), "only the changes in the coefficients are tested")
})

test_that("the sequential statistic is NA when no segment can hold another break", {
  fit3 <- breaks_sur(y1, sur_x, m = 3, trim = 0.2)
  expect_warning(tb <- test_breaks(fit3), "is NA")
  expect_true(is.na(tb$seq["Seq(3|2)", "statistic"]))
  expect_false(is.na(tb$seq["Seq(2|1)", "statistic"]))
})

test_that("select_breaks chooses two breaks on the example and refits", {
  fit3 <- breaks_sur(y1, sur_x, m = 3, trim = 0.2)
  sel <- suppressWarnings(select_breaks(fit3))
  expect_s3_class(sel, "mvbreakSelect")
  expect_equal(sel$n_breaks, 2L)
  expect_equal(sel$fit$dates, c(56L, 100L))
  expect_equal(sel$first_test, "WDmax")
  sel2 <- suppressWarnings(select_breaks(fit3, first_test = "supLR"))
  expect_equal(sel2$n_breaks, 2L)
  tab <- summary(sel)
  expect_s3_class(tab, "data.frame")
  expect_true(all(c("statistic", "reject_at") %in% names(tab)))
})

test_that("tests argument and max_breaks are validated", {
  expect_error(test_breaks(fit2, tests = "foo"), "'tests' must be")
  expect_error(test_breaks(fit2, max_breaks = 0), "'max_breaks' must be a positive integer")
  tb <- test_breaks(fit2, tests = "supLR")
  expect_null(tb$WDmax); expect_null(tb$seq)
  expect_output(print(tb), "supLR test of no break")
})

test_that("a user-supplied R allows the test at the fitted m only", {
  S <- rbind(c(1,0,0,0), c(0,1,0,0), c(0,0,0,0), c(0,0,1,0), c(0,0,0,0), c(0,0,0,1))
  R <- mvbreak:::restriction_matrix(4, 2, c(1, 3))
  fit <- breaks_system(y1, z1, S, m = 2, trim = 0.2, R = R)
  expect_true(fit$settings$R_user)
  expect_message(tb <- test_breaks(fit), "only the supLR test for 2 break")
  expect_equal(unname(tb$supLR[1, "statistic"]), 48.382, tolerance = 1e-4)
  expect_error(test_breaks(fit, max_breaks = 3), "must be 2")
  expect_error(select_breaks(fit), "needs restrictions given through 'fixed'")
})
