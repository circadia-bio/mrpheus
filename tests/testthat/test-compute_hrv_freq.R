# tests/testthat/test-compute_hrv_freq.R
#
# Fully synthetic tests for compute_hrv_freq() (new in mrpheus 0.1.4.9000).
# No participant-derived data is used, per project convention.

.make_synthetic_qrs <- function(dur, base_rr, fs, resp_freq = NULL, rsa_amp = 0) {
  t     <- 0
  times <- numeric(0)
  while (t < dur) {
    rr <- base_rr
    if (!is.null(resp_freq)) rr <- rr + rsa_amp * sin(2 * pi * resp_freq * t)
    t <- t + rr
    if (t < dur) times <- c(times, t)
  }
  structure(
    list(
      qrs_i   = as.integer(round(times * fs)),
      qrs_amp = rep(1, length(times)),
      delay   = 0,
      fs      = fs
    ),
    class = "mrpheus_qrs"
  )
}

.make_resp_signal <- function(dur, fs, freq, noise_sd = 0.1, seed = 1) {
  set.seed(seed)
  t <- seq(0, dur, by = 1 / fs)
  sin(2 * pi * freq * t) + rnorm(length(t), sd = noise_sd)
}

test_that("compute_hrv_freq recovers a known injected respiratory frequency", {
  fs  <- 250
  dur <- 180
  f0  <- 0.9

  qrs  <- .make_synthetic_qrs(dur, base_rr = 60 / 140, fs = fs,
                               resp_freq = f0, rsa_amp = 0.04)
  resp <- .make_resp_signal(dur, fs, f0, seed = 1)

  out <- compute_hrv_freq(qrs, resp, fs)

  expect_equal(out$resp_freq_hz, f0, tolerance = 0.05)
  expect_true(is.finite(out$hf_power))
  expect_gt(out$hf_power, 0)
})

test_that("compute_hrv_freq recovers a different injected frequency", {
  fs  <- 250
  dur <- 180
  f0  <- 1.4  # faster neonatal breathing rate (~84 breaths/min)

  qrs  <- .make_synthetic_qrs(dur, base_rr = 60 / 150, fs = fs,
                               resp_freq = f0, rsa_amp = 0.03)
  resp <- .make_resp_signal(dur, fs, f0, seed = 2)

  out <- compute_hrv_freq(qrs, resp, fs)

  expect_equal(out$resp_freq_hz, f0, tolerance = 0.05)
})

test_that("compute_hrv_freq: HF power increases with RSA modulation amplitude", {
  fs  <- 250
  dur <- 180
  f0  <- 0.9

  resp <- .make_resp_signal(dur, fs, f0, seed = 3)

  qrs_low  <- .make_synthetic_qrs(dur, base_rr = 60 / 140, fs = fs,
                                   resp_freq = f0, rsa_amp = 0.01)
  qrs_high <- .make_synthetic_qrs(dur, base_rr = 60 / 140, fs = fs,
                                   resp_freq = f0, rsa_amp = 0.08)

  out_low  <- compute_hrv_freq(qrs_low,  resp, fs)
  out_high <- compute_hrv_freq(qrs_high, resp, fs)

  expect_gt(out_high$hf_power, out_low$hf_power)
})

test_that("compute_hrv_freq returns band bounds centered on the detected frequency", {
  fs  <- 250
  dur <- 180
  f0  <- 0.9

  qrs  <- .make_synthetic_qrs(dur, base_rr = 60 / 140, fs = fs,
                               resp_freq = f0, rsa_amp = 0.04)
  resp <- .make_resp_signal(dur, fs, f0, seed = 4)

  out <- compute_hrv_freq(qrs, resp, fs, band_halfwidth_hz = 0.15)

  expect_equal(out$band_high - out$band_low, 0.3, tolerance = 1e-8)
  expect_equal((out$band_low + out$band_high) / 2, out$resp_freq_hz, tolerance = 1e-8)
})

test_that("compute_hrv_freq returns NA fields when there are too few beats", {
  qrs <- structure(
    list(qrs_i = as.integer(c(100, 200, 300)), qrs_amp = c(1, 1, 1),
         delay = 0, fs = 250),
    class = "mrpheus_qrs"
  )
  resp <- .make_resp_signal(10, 250, 0.9, seed = 5)

  out <- compute_hrv_freq(qrs, resp, fs = 250)

  expect_true(is.na(out$hf_power))
  expect_true(is.na(out$resp_freq_hz))
  expect_equal(out$n_beats_used, 2L)  # 3 peaks -> 2 RR intervals, both < 10 -> empty result
})

test_that("compute_hrv_freq drops physiologically implausible RR intervals", {
  fs <- 250
  qrs_times <- c(0.4, 0.8, 1.2, 5.0, 5.4, 5.8, 6.2, 6.6, 7.0, 7.4, 7.8, 8.2)
  # the jump from 1.2s to 5.0s is a ~3.8s RR interval (outside 60-220 bpm)
  qrs <- structure(
    list(
      qrs_i   = as.integer(round(qrs_times * fs)),
      qrs_amp = rep(1, length(qrs_times)),
      delay   = 0,
      fs      = fs
    ),
    class = "mrpheus_qrs"
  )
  resp <- .make_resp_signal(9, fs, 0.9, seed = 6)

  out <- compute_hrv_freq(qrs, resp, fs)
  # 12 peaks -> 11 raw RR intervals, one implausible -> at most 10 valid
  expect_lte(out$n_beats_used, 10L)
})
