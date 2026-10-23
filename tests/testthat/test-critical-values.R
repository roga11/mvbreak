# Critical values: the shipped table against Bai and Perron's published values, the
# consistency between table and simulation, and the seed handling.

test_that("the table reproduces Bai and Perron (1998) for one changing parameter", {
  cv <- critical_values(q = 1, m = 5, trm = 0.05)
  expect_equal(rownames(cv), c("90%", "95%", "97.5%", "99%"))
  expect_equal(cv["95%", "supLR(m=1)"], 9.63, tolerance = 0.03)      # Table I, sup F(1), q = 1
  expect_equal(cv["95%", "supLR(m=2)"] / 2, 8.78, tolerance = 0.03)  # their tables report supLR(m)/m
  expect_equal(cv["95%", "UDmax"], 10.17, tolerance = 0.03)          # double maximum tests with M = 5
  expect_equal(cv["95%", "WDmax"], 10.91, tolerance = 0.03)
  expect_equal(attr(cv, "method"), "table")
})

test_that("values increase with q and with the number of breaks, and decrease with trimming", {
  a <- critical_values(q = 2, m = 3, trm = 0.15)
  b <- critical_values(q = 6, m = 3, trm = 0.15)
  expect_true(all(b[, 1:3] > a[, 1:3]))
  expect_true(all(diff(a["95%", paste0("supLR(m=", 1:3, ")")]) > 0))
  c15 <- critical_values(q = 2, m = 1, trm = 0.15)["95%", "supLR(m=1)"]
  c25 <- critical_values(q = 2, m = 1, trm = 0.25)["95%", "supLR(m=1)"]
  expect_true(c25 < c15)
})

test_that("a simulated value is close to the table and a seed restores the RNG state", {
  skip_on_cran()
  set.seed(42); before <- .Random.seed
  cv <- critical_values(q = 2, m = 1, trm = 0.15, method = "simulate", B = 2000, n = 300, seed = 7)
  expect_identical(.Random.seed, before)
  expect_equal(attr(cv, "method"), "simulate")
  tab <- critical_values(q = 2, m = 1, trm = 0.15)
  expect_equal(cv["95%", "supLR(m=1)"], tab["95%", "supLR(m=1)"], tolerance = 0.1)
  # same seed, same draws
  cv2 <- critical_values(q = 2, m = 1, trm = 0.15, method = "simulate", B = 2000, n = 300, seed = 7)
  expect_identical(cv, cv2)
})

test_that("cases outside the table are simulated with a message", {
  skip_on_cran()
  expect_message(cv <- critical_values(q = 2, m = 1, trm = 0.12, B = 500, n = 200, seed = 1), "does not cover")
  expect_equal(attr(cv, "method"), "simulate")
})

test_that("inputs are validated", {
  expect_error(critical_values(q = 0), "'q' must be a positive integer")
  expect_error(critical_values(q = 2, m = 4, trm = 0.25), "at most 2 break")
  expect_error(critical_values(q = 2, level = 1), "'level' must lie")
})
