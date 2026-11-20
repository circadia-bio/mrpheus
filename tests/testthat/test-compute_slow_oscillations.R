# tests/testthat/test-compute_slow_oscillations.R
#
# Behaviour of compute_slow_oscillations() and its helpers. Agreement with
# yasa.sw_detect() is tested in test-yasa-slow-oscillation-parity.R.
#
# All fixtures are fully synthetic (no participant-derived data).
#
# The oscillation is -sin on purpose: the filtered epoch starts negative, which
# used to make the old zero-crossing detector pair crossings half a cycle out of
# phase and return no events.

.make_so_psg <- function(n_epochs = 3L, sr = 256, freq = 1, amp_uv = 100,
                         epoch_s = 30) {
  n  <- as.integer(epoch_s * sr)
  tt <- (seq_len(n) - 1L) / sr
  structure(
    list(
      epochs      = replicate(n_epochs, list(C3 = -sin(2 * pi * freq * tt) * amp_uv),
                              simplify = FALSE),
      n_epochs    = n_epochs,
      epoch_s     = epoch_s,
      channel_map = tibble::tibble(label = "C3", type = "EEG",
                                   sample_rate = sr, bad = FALSE)
    ),
    class = "mrpheus_psg"
  )
}

test_that("compute_slow_oscillations runs at 256 Hz and finds 1 Hz oscillations", {
  so <- compute_slow_oscillations(.make_so_psg(sr = 256), channel = "C3")
  expect_s3_class(so, "tbl_df")
  expect_gt(nrow(so), 0L)
  # defaults are YASA's: ptp 75-350, trough 40-200, peak 10-150 uV
  expect_true(all(so$ptp_uv >= 75 & so$ptp_uv <= 350))
  expect_true(all(-so$neg_peak_uv >= 40 & -so$neg_peak_uv <= 200))
  expect_true(all(so$pos_peak_uv >= 10 & so$pos_peak_uv <= 150))
})

test_that("compute_slow_oscillations runs at 100 Hz too", {
  so <- compute_slow_oscillations(.make_so_psg(sr = 100), channel = "C3")
  expect_gt(nrow(so), 0L)
})

test_that("trough and peak limits are applied, and NULL switches them off", {
  psg <- .make_so_psg(sr = 256)           # about -/+ 97 uV after filtering
  expect_equal(nrow(compute_slow_oscillations(psg, amp_neg_uv = c(150, 200))), 0L)
  expect_equal(nrow(compute_slow_oscillations(psg, amp_pos_uv = c(120, 150))), 0L)
  expect_gt(nrow(compute_slow_oscillations(psg, amp_neg_uv = NULL, amp_pos_uv = NULL)), 0L)
})

test_that("a transition band of zero is rejected", {
  expect_error(compute_slow_oscillations(.make_so_psg(sr = 256), transition_hz = 0),
               "Invalid slow oscillation filter band")
})

test_that("stages restricts the epochs in which peaks are accepted", {
  psg <- .make_so_psg(n_epochs = 4L, sr = 100)
  all <- compute_slow_oscillations(psg)
  one <- compute_slow_oscillations(psg, stages = 2L)
  expect_gt(nrow(all), nrow(one))
  expect_gt(nrow(one), 0L)
  expect_true(all(one$epoch == 2L))
  expect_equal(one$neg_peak_s, all$neg_peak_s[all$epoch == 2L])
})

test_that("epoch and the times relative to it are consistent", {
  so <- compute_slow_oscillations(.make_so_psg(n_epochs = 3L, sr = 256))
  expect_equal(so$epoch, as.integer(floor(so$neg_peak_s / 30)) + 1L)
  expect_equal(so$duration_s, so$end_s - so$start_s, tolerance = 1e-3)
  expect_true(all(so$start_s < (so$neg_peak_s - (so$epoch - 1) * 30)))
  expect_true(all(so$neg_peak_s < so$mid_crossing_s & so$mid_crossing_s < so$pos_peak_s))
  expect_true(all(so$slope_uv_s > 0))
})

test_that("a channel in the wrong unit is skipped with a warning", {
  psg <- .make_so_psg(sr = 100)
  psg$epochs <- lapply(psg$epochs, function(e) list(C3 = e$C3 * 1e-6))   # volts, not uV
  expect_warning(so <- compute_slow_oscillations(psg), "Wrong data amplitude")
  expect_equal(nrow(so), 0L)
  expect_named(so, c("epoch", "start_s", "end_s", "duration_s", "neg_peak_uv", "pos_peak_uv",
                     "ptp_uv", "channel", "neg_peak_s", "mid_crossing_s", "pos_peak_s", "slope_uv_s"))
})

test_that("a quiet channel returns an empty result without a warning", {
  set.seed(1)
  psg <- .make_so_psg(sr = 100)
  psg$epochs <- lapply(psg$epochs, function(e) list(C3 = rnorm(3000)))
  expect_no_warning(so <- suppressMessages(compute_slow_oscillations(psg)))
  expect_equal(nrow(so), 0L)
})

test_that(".sw_find_peaks follows scipy.signal.find_peaks", {
  x <- c(0, 1, 0, 2, 2, 0, 3, 3, 3, 1)         # scipy: peaks 1, 3 (plateau of 2), 7 (plateau of 3)
  expect_equal(mrpheus:::.sw_find_peaks(x), c(1L, 3L, 7L))
  expect_equal(mrpheus:::.sw_find_peaks(x, c(1.5, 2.5)), 3L)              # height limits, inclusive
  expect_equal(mrpheus:::.sw_find_peaks(x, c(2, 3)), c(3L, 7L))
  expect_equal(mrpheus:::.sw_find_peaks(c(2, 1, 0, 1, 2)), integer(0))   # edges are never peaks
  expect_equal(mrpheus:::.sw_find_peaks(c(0, 1, 1)), integer(0))         # a plateau that never drops
  expect_equal(mrpheus:::.sw_find_peaks(c(0, 1)), integer(0))
})

test_that("the MNE FIR design has the requested band and length", {
  for (fs in c(100, 256)) {
    h <- mrpheus:::.mne_fir_bandpass_design(fs, 0.3, 1.5, 0.2, 0.2)
    expect_equal(length(h), as.integer(round(3.3 / 0.2 * fs)) + 1L)         # odd: linear phase
    gain <- function(f) Mod(sum(h * exp(-2i * pi * f * (seq_along(h) - 1) / fs)))
    expect_lt(gain(0), 1e-3)                                                # no DC
    expect_equal(gain(0.9), 1, tolerance = 1e-3)                            # pass band
    expect_equal(gain(0.2), 0.5, tolerance = 0.02)                          # -6 dB at the cutoffs
    expect_equal(gain(1.6), 0.5, tolerance = 0.02)
    expect_lt(gain(3), 1e-3)                                                # stop band
  }
})

test_that("an impossible filter band is rejected", {
  psg <- .make_so_psg(sr = 256)
  # 0.05 - 0.2 < 0: the lower stop frequency of the FIR filter would be negative
  expect_error(compute_slow_oscillations(psg, freq_so = c(0.05, 1.5)),
               "Invalid slow oscillation filter band")
  expect_error(compute_slow_oscillations(psg, freq_so = c(0.3, 200)),
               "Invalid slow oscillation filter band")
})
