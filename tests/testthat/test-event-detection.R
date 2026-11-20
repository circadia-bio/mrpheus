# tests/testthat/test-event-detection.R
#
# Unit + parity tests for the Rcpp functions of the staging normalisation
# pipeline: roll_right_mean_cpp and robust_scale_cpp, and rowmedian_cpp from
# .welch_median_psd. (The event detection functions roll_rms_cpp and
# detect_so_candidates_cpp were removed when compute_spindles() and
# compute_slow_oscillations() were ported from YASA; see their parity tests.)

# ── Pure-R references (inline — no zoo at test time) ─────────────────────────

.roll_right_r <- function(x, k) {
  n <- length(x)
  vapply(seq_len(n), function(i) {
    mean(x[max(1L, i - k + 1L):i])
  }, numeric(1))
}

.robust_r <- function(x, q_low = 0.05, q_high = 0.95) {
  xs  <- sort(x[!is.na(x)])
  n   <- length(xs)
  med <- if (n %% 2 == 1) xs[n %/% 2 + 1] else (xs[n %/% 2] + xs[n %/% 2 + 1]) / 2
  q7  <- function(p) {
    h <- (n - 1) * p; lo <- floor(h) + 1L; hi <- lo + 1L
    if (hi > n) return(xs[n])
    xs[lo] + (h - floor(h)) * (xs[hi] - xs[lo])
  }
  scale <- q7(q_high) - q7(q_low) + 1e-10
  ifelse(is.na(x), NA_real_, (x - med) / scale)
}


# ── rowmedian_cpp ─────────────────────────────────────────────────────────────

test_that("rowmedian_cpp returns correct length", {
  m <- matrix(rnorm(251 * 9), nrow = 251)
  expect_length(rowmedian_cpp(m), 251L)
})

test_that("rowmedian_cpp matches apply(m, 1, median) on odd ncol", {
  set.seed(1)
  m <- matrix(rnorm(251 * 9), nrow = 251)
  expect_equal(rowmedian_cpp(m), apply(m, 1L, median), tolerance = 1e-12)
})

test_that("rowmedian_cpp matches apply(m, 1, median) on even ncol", {
  set.seed(2)
  m <- matrix(rnorm(100 * 8), nrow = 100)
  expect_equal(rowmedian_cpp(m), apply(m, 1L, median), tolerance = 1e-12)
})

test_that("rowmedian_cpp returns constant row for constant matrix", {
  m   <- matrix(rep(5, 20), nrow = 4)
  out <- rowmedian_cpp(m)
  expect_equal(out, rep(5, 4), tolerance = 1e-12)
})

# ── roll_right_mean_cpp ───────────────────────────────────────────────────────

test_that("roll_right_mean_cpp length matches input", {
  x <- rnorm(50)
  expect_length(roll_right_mean_cpp(x, 4L), 50L)
})

test_that("roll_right_mean_cpp k=1 is identity", {
  set.seed(5)
  x <- rnorm(30)
  expect_equal(roll_right_mean_cpp(x, 1L), x, tolerance = 1e-12)
})

test_that("roll_right_mean_cpp partial windows at start are correct", {
  x <- c(2, 4, 6, 8, 10)
  out <- roll_right_mean_cpp(x, 4L)
  expect_equal(out[1], 2,              tolerance = 1e-12)  # window [1]
  expect_equal(out[2], mean(x[1:2]),   tolerance = 1e-12)  # window [1:2]
  expect_equal(out[3], mean(x[1:3]),   tolerance = 1e-12)  # window [1:3]
  expect_equal(out[4], mean(x[1:4]),   tolerance = 1e-12)  # full window
  expect_equal(out[5], mean(x[2:5]),   tolerance = 1e-12)  # sliding
})

test_that("roll_right_mean_cpp matches R reference on random signal", {
  set.seed(6)
  x <- rnorm(200)
  expect_equal(roll_right_mean_cpp(x, 4L), .roll_right_r(x, 4L), tolerance = 1e-12)
})

# ── robust_scale_cpp ──────────────────────────────────────────────────────────

test_that("robust_scale_cpp length matches input", {
  expect_length(robust_scale_cpp(rnorm(100)), 100L)
})

test_that("robust_scale_cpp median of output is 0", {
  set.seed(7)
  x <- rnorm(500)
  expect_equal(median(robust_scale_cpp(x)), 0, tolerance = 1e-6)
})

test_that("robust_scale_cpp constant input returns zeros", {
  expect_equal(unique(robust_scale_cpp(rep(5, 50))), 0, tolerance = 1e-5)
})

test_that("robust_scale_cpp preserves NA positions", {
  x   <- c(1, 2, NA, 4, 5)
  out <- robust_scale_cpp(x)
  expect_true(is.na(out[3]))
  expect_false(anyNA(out[-3]))
})

test_that("robust_scale_cpp matches R reference on random signal", {
  set.seed(8)
  x <- rnorm(1000)
  expect_equal(robust_scale_cpp(x), .robust_r(x), tolerance = 1e-10)
})

test_that("robust_scale_cpp respects custom q_low and q_high", {
  set.seed(9)
  x <- rnorm(500)
  expect_equal(robust_scale_cpp(x, 0.1, 0.9), .robust_r(x, 0.1, 0.9), tolerance = 1e-10)
})
