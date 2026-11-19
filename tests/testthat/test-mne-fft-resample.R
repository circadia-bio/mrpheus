# tests/testthat/test-mne-fft-resample.R
#
# Tests for `.mne_fft_resample()`, the internal FFT resampler used by
# stage_epochs() to bring EEG/EOG/EMG to 100 Hz before feature extraction.
#
# The downsampling branch (new_len < padded length) was missing: staging a
# 256 Hz recording failed with "only 0's may be mixed with negative
# subscripts". The expected values below were checked against
# mne.filter.resample(x, up = 100/256, down = 1, npad = "auto", method = "fft")
# in Python (MNE 1.13.2): maximum absolute difference ~5e-14 on four signal
# lengths. The test here uses an analytic sinusoid so that it needs no Python.
#
# All fixtures are fully synthetic (no participant-derived data).

test_that(".mne_fft_resample downsamples 256 Hz to 100 Hz with the right length", {
  for (n in c(7680L, 123457L, 100001L)) {
    x <- sin(2 * pi * 10 * (seq_len(n) - 1L) / 256)
    y <- .mne_fft_resample(x, up = 100 / 256)
    expect_length(y, as.integer(round(n * 100 / 256)))
  }
})

test_that(".mne_fft_resample recovers a 10 Hz sinusoid when downsampling 256 -> 100 Hz", {
  n  <- 7680L                                   # 30 s, 300 whole cycles
  x  <- sin(2 * pi * 10 * (seq_len(n) - 1L) / 256)
  y  <- .mne_fft_resample(x, up = 100 / 256)
  ref <- sin(2 * pi * 10 * (seq_along(y) - 1L) / 100)

  # Observed in the reference check: ~9e-4 at the edges, ~3e-6 in the interior.
  expect_lt(max(abs(y - ref)), 5e-3)
  interior <- 101:(length(y) - 100)
  expect_lt(max(abs(y[interior] - ref[interior])), 1e-4)
})

test_that(".mne_fft_resample downsampling preserves a DC offset", {
  n <- 7680L
  x <- 50 + sin(2 * pi * 10 * (seq_len(n) - 1L) / 256)
  y <- .mne_fft_resample(x, up = 100 / 256)
  expect_equal(mean(y[101:(length(y) - 100)]), 50, tolerance = 1e-3)
})

test_that(".mne_fft_resample upsampling path is unchanged", {
  n <- 6000L
  x <- sin(2 * pi * 5 * (seq_len(n) - 1L) / 100)
  y <- .mne_fft_resample(x, up = 2.5)
  expect_length(y, as.integer(round(n * 2.5)))
  ref <- sin(2 * pi * 5 * (seq_along(y) - 1L) / 250)
  interior <- 251:(length(y) - 250)
  expect_lt(max(abs(y[interior] - ref[interior])), 1e-3)
})
