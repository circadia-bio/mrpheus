# tests/testthat/test-compute_slow_oscillations.R
#
# Regression test for compute_slow_oscillations() at 256 Hz. The 0.5-2 Hz
# band-pass was designed as one 8th-order transfer function, which is
# numerically singular at 256 Hz (gsignal::filter_zi: rcond ~ 4e-18), so every
# call failed with "system is computationally singular". It passed at lower
# sample rates, and the existing tests only exercise detect_so_candidates_cpp(),
# not the filtering step.
#
# All fixtures are fully synthetic (no participant-derived data).
#
# The oscillation is -sin on purpose: after filtering, the epoch starts
# negative, which used to make detect_so_candidates_cpp() pair zero crossings
# half a cycle out of phase and return no events (see test-event-detection.R).

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

test_that("the old, wider criteria are still available through the arguments", {
  psg <- .make_so_psg(sr = 256)
  so  <- compute_slow_oscillations(
    psg, freq_so = c(0.5, 2), amp_ptp_threshold_uv = c(75, 500),
    duration_neg_s = c(0.1, 1.5), amp_neg_uv = NULL, amp_pos_uv = NULL,
    transition_hz = 0
  )
  expect_gt(nrow(so), 0L)
})

test_that("an impossible filter band is rejected", {
  psg <- .make_so_psg(sr = 256)
  # 0.05 - 0.2 / 2 < 0: the high-pass cutoff would be negative
  expect_error(compute_slow_oscillations(psg, freq_so = c(0.05, 1.5)),
               "Invalid slow oscillation filter band")
  expect_error(compute_slow_oscillations(psg, freq_so = c(0.3, 200)),
               "Invalid slow oscillation filter band")
})
