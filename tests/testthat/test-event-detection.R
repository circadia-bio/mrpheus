# tests/testthat/test-event-detection.R
#
# Unit + parity tests for event detection Rcpp functions:
#   roll_rms_cpp              (compute_spindles hot path)
#   detect_so_candidates_cpp  (compute_slow_oscillations hot path)
#
# Also covers roll_right_mean_cpp and robust_scale_cpp from the staging
# normalisation pipeline, and rowmedian_cpp from .welch_median_psd.

# ── Pure-R references (inline — no zoo at test time) ─────────────────────────

.roll_rms_r <- function(x, k) {
  n     <- length(x)
  left  <- (k - 1L) %/% 2L
  right <- k - 1L - left
  out   <- rep(NA_real_, n)
  for (i in seq(left + 1L, n - right)) {
    w <- x[(i - left):(i + right)]
    out[i] <- sqrt(mean(w^2))
  }
  out
}

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

.so_r <- function(sig, sr, dur_neg, dur_pos, amp) {
  zc <- which(diff(sign(sig)) != 0)
  if (length(zc) < 4) return(matrix(numeric(0), ncol = 6))
  rows <- lapply(seq(1, length(zc) - 3, by = 2), function(k) {
    ns <- zc[k]; ne <- zc[k + 1]; pe <- zc[k + 2]
    dn <- (ne - ns) / sr; dp <- (pe - ne) / sr
    if (dn < dur_neg[1] || dn > dur_neg[2]) return(NULL)
    if (dp < dur_pos[1] || dp > dur_pos[2]) return(NULL)
    np <- min(sig[ns:ne]); pp <- max(sig[ne:pe]); ptp <- pp - np
    if (ptp < amp[1] || ptp > amp[2]) return(NULL)
    c(ns, ne, pe, np, pp, ptp)
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(matrix(numeric(0), ncol = 6))
  do.call(rbind, rows)
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

# ── roll_rms_cpp ──────────────────────────────────────────────────────────────

test_that("roll_rms_cpp length matches input", {
  expect_length(roll_rms_cpp(rnorm(100), 5L), 100L)
})

test_that("roll_rms_cpp edge positions are NA", {
  out <- roll_rms_cpp(rnorm(50), 5L)
  # left = floor((5-1)/2) = 2; right = 2; positions 1,2 and 49,50 are NA
  expect_true(all(is.na(out[1:2])))
  expect_true(all(is.na(out[49:50])))
})

test_that("roll_rms_cpp is non-negative for all non-NA values", {
  out <- roll_rms_cpp(rnorm(200), 15L)
  expect_true(all(out[!is.na(out)] >= 0))
})

test_that("roll_rms_cpp constant input: rms = |constant|", {
  out <- roll_rms_cpp(rep(-3, 100), 5L)
  expect_equal(out[3:98], rep(3, 96), tolerance = 1e-10)
})

test_that("roll_rms_cpp matches R reference for odd k", {
  set.seed(10)
  x <- rnorm(300)
  expect_equal(roll_rms_cpp(x, 15L), .roll_rms_r(x, 15L), tolerance = 1e-10)
})

test_that("roll_rms_cpp matches R reference for even k", {
  set.seed(11)
  x <- rnorm(300)
  expect_equal(roll_rms_cpp(x, 30L), .roll_rms_r(x, 30L), tolerance = 1e-10)
})

# ── detect_so_candidates_cpp ──────────────────────────────────────────────────

test_that("detect_so_candidates_cpp returns 0 rows for flat signal", {
  out <- detect_so_candidates_cpp(rep(0, 500), 256, 0.1, 1.5, 0.1, 1.0, 50, 500)
  expect_equal(nrow(out), 0L)
})

test_that("detect_so_candidates_cpp returns 0 rows for < 4 zero crossings", {
  out <- detect_so_candidates_cpp(c(1, -1, 1), 256, 0.1, 1.5, 0.1, 1.0, 0, 1e6)
  expect_equal(nrow(out), 0L)
})

test_that("detect_so_candidates_cpp output has 6 columns with correct names", {
  sr  <- 256
  sig <- -sin(2 * pi * 1 * seq(0, 4, by = 1/sr)) * 200
  out <- detect_so_candidates_cpp(sig, sr, 0.1, 1.5, 0.1, 1.0, 50, 500)
  expect_equal(ncol(out), 6L)
  expect_equal(colnames(out),
               c("neg_start", "neg_end", "pos_end", "neg_peak", "pos_peak", "ptp"))
})

test_that("detect_so_candidates_cpp ptp = pos_peak - neg_peak", {
  sr  <- 256
  sig <- -sin(2 * pi * 1 * seq(0, 4, by = 1/sr)) * 200
  out <- detect_so_candidates_cpp(sig, sr, 0.1, 1.5, 0.1, 1.0, 50, 500)
  if (nrow(out) > 0) {
    expect_equal(out[, "ptp"], out[, "pos_peak"] - out[, "neg_peak"], tolerance = 1e-10)
  }
})

test_that("detect_so_candidates_cpp ptp is always within amp bounds", {
  sr  <- 256
  sig <- -sin(2 * pi * 1 * seq(0, 4, by = 1/sr)) * 200
  out <- detect_so_candidates_cpp(sig, sr, 0.1, 1.5, 0.1, 1.0, 50, 500)
  if (nrow(out) > 0) {
    expect_true(all(out[, "ptp"] >= 50 & out[, "ptp"] <= 500))
  }
})

test_that("detect_so_candidates_cpp indices are 1-indexed and ordered", {
  sr  <- 256
  sig <- -sin(2 * pi * 1 * seq(0, 4, by = 1/sr)) * 200
  out <- detect_so_candidates_cpp(sig, sr, 0.1, 1.5, 0.1, 1.0, 50, 500)
  if (nrow(out) > 0) {
    expect_true(all(out[, 1] >= 1))
    expect_true(all(out[, 1] < out[, 2]))
    expect_true(all(out[, 2] < out[, 3]))
  }
})

test_that("detect_so_candidates_cpp matches R reference (count and indices)", {
  set.seed(42)
  sr  <- 256
  sig <- -sin(2 * pi * 1 * seq(0, 8, by = 1/sr)) * 200 + rnorm(256 * 8 + 1, sd = 10)
  cpp <- detect_so_candidates_cpp(sig, sr, 0.1, 1.5, 0.1, 1.0, 50, 500)
  r   <- .so_r(sig, sr, c(0.1, 1.5), c(0.1, 1.0), c(50, 500))
  expect_equal(nrow(cpp), nrow(r))
  if (nrow(cpp) > 0 && nrow(r) > 0) {
    expect_equal(cpp[, 1], r[, 1], tolerance = 1e-10)  # neg_start
    expect_equal(cpp[, 3], r[, 3], tolerance = 1e-10)  # pos_end
    expect_equal(cpp[, 6], r[, 6], tolerance = 1e-10)  # ptp
  }
})

test_that("detect_so_candidates_cpp amp threshold filters correctly", {
  sr  <- 256
  # Large amplitude signal: all candidates should pass loose bounds, fail tight ones
  sig <- -sin(2 * pi * 1 * seq(0, 4, by = 1/sr)) * 500
  loose <- detect_so_candidates_cpp(sig, sr, 0.1, 1.5, 0.1, 1.0, 1, 2000)
  tight <- detect_so_candidates_cpp(sig, sr, 0.1, 1.5, 0.1, 1.0, 1, 10)
  expect_gte(nrow(loose), nrow(tight))
})
