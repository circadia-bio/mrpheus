# tests/testthat/test-boxcar-integration.R
#
# Equivalence test for the cumulative-sum boxcar moving-average used in
# detect_qrs() Stage 4, replacing stats::convolve() (see NEWS.md, mrpheus
# 0.1.4.9000). The boxcar step is inlined inside detect_qrs() rather than a
# standalone function, so this test reproduces the exact formula used there
# and checks it against the FFT-based stats::convolve() reference it
# replaced, on fully synthetic vectors of varying length. A final smoke test
# runs the full detect_qrs() pipeline end to end.

.cumsum_boxcar <- function(x, win) {
  n_s     <- length(x)
  padded  <- c(rep(0, win - 1), x, rep(0, win - 1))
  cs      <- cumsum(c(0, padded))
  out_len <- n_s + win - 1
  (cs[(win + 1):(win + out_len)] - cs[1:out_len]) / win
}

test_that("cumsum boxcar matches stats::convolve(type = 'open') on odd-length signal", {
  set.seed(201)
  x   <- rnorm(1001)
  win <- 31L
  ref <- as.double(stats::convolve(x, rep(1, win) / win, type = "open"))
  out <- .cumsum_boxcar(x, win)
  expect_equal(out, ref, tolerance = 1e-10)
})

test_that("cumsum boxcar matches stats::convolve(type = 'open') on even-length signal", {
  set.seed(202)
  x   <- rnorm(2000)
  win <- 40L
  ref <- as.double(stats::convolve(x, rep(1, win) / win, type = "open"))
  out <- .cumsum_boxcar(x, win)
  expect_equal(out, ref, tolerance = 1e-10)
})

test_that("cumsum boxcar matches stats::convolve for a range of window sizes", {
  set.seed(203)
  x <- rnorm(500)
  for (win in c(3L, 7L, 15L, 30L, 75L)) {
    ref <- as.double(stats::convolve(x, rep(1, win) / win, type = "open"))
    out <- .cumsum_boxcar(x, win)
    expect_equal(out, ref, tolerance = 1e-10,
                 info = sprintf("mismatch at win = %d", win))
  }
})

test_that("cumsum boxcar output length matches full convolution length", {
  x   <- rnorm(150)
  win <- 25L
  out <- .cumsum_boxcar(x, win)
  expect_equal(length(out), length(x) + win - 1L)
})

test_that("detect_qrs() completes quickly on a longer synthetic signal", {
  # Regression test for the *combined* effect of both Stage 4/5 fixes: a
  # synthetic signal long enough (~2 min) that the pre-fix FFT convolution
  # or O(n*k) peak selection would have been noticeably slow or hanging.
  set.seed(204)
  fs      <- 200
  dur     <- 120  # seconds
  hr_hz   <- 2    # 120 bpm
  t       <- seq(0, dur, by = 1 / fs)
  rr      <- 1 / hr_hz
  beat_t  <- seq(0.5, dur - 0.5, by = rr)

  pulse_width <- 0.02
  ecg <- numeric(length(t))
  for (pt in beat_t) {
    ecg <- ecg + exp(-((t - pt)^2) / (2 * pulse_width^2))
  }
  ecg <- ecg + rnorm(length(t), sd = 0.01)

  t0  <- Sys.time()
  qrs <- detect_qrs(ecg, fs = fs)
  elapsed <- as.numeric(Sys.time() - t0, units = "secs")

  expect_s3_class(qrs, "mrpheus_qrs")
  expect_lt(elapsed, 15)
  expect_gt(length(qrs$qrs_i), 0L)
})
