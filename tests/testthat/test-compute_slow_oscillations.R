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
# NOTE on the sign of the test oscillation: detect_so_candidates_cpp() pairs
# zero crossings in steps of two and assumes the first crossing in the epoch
# starts a NEGATIVE half-wave, without checking polarity. A filtered epoch that
# starts negative (e.g. -sin) therefore yields zero events, while +sin yields
# events. The existing C++ tests hide this because their signals start at
# exactly 0. The oscillation below is +sin so that this file tests the
# filtering fix only; the polarity issue is a separate, open problem.

.make_so_psg <- function(n_epochs = 3L, sr = 256, freq = 1, amp_uv = 100,
                         epoch_s = 30) {
  n  <- as.integer(epoch_s * sr)
  tt <- (seq_len(n) - 1L) / sr
  structure(
    list(
      epochs      = replicate(n_epochs, list(C3 = sin(2 * pi * freq * tt) * amp_uv),
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
  expect_true(all(so$ptp_uv >= 75 & so$ptp_uv <= 500))
})

test_that("compute_slow_oscillations runs at 100 Hz too", {
  so <- compute_slow_oscillations(.make_so_psg(sr = 100), channel = "C3")
  expect_gt(nrow(so), 0L)
})
