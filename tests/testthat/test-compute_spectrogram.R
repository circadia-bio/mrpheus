# tests/testthat/test-compute_spectrogram.R
#
# Regression test for compute_spectrogram(). It called gsignal::pwelch() with
# `noverlap`, which gsignal does not have (it takes `overlap` as a fraction),
# so every call failed with "unused argument (noverlap = novlp)".
#
# All fixtures are fully synthetic (no participant-derived data).

.make_sine_psg <- function(n_epochs = 3L, sr = 100, freq = 10, epoch_s = 30) {
  n  <- as.integer(epoch_s * sr)
  tt <- (seq_len(n) - 1L) / sr
  structure(
    list(
      epochs      = replicate(n_epochs, list(C3 = sin(2 * pi * freq * tt)),
                              simplify = FALSE),
      n_epochs    = n_epochs,
      epoch_s     = epoch_s,
      channel_map = tibble::tibble(label = "C3", type = "EEG",
                                   sample_rate = sr, bad = FALSE)
    ),
    class = "mrpheus_psg"
  )
}

test_that("compute_spectrogram runs and returns one row per epoch", {
  sg <- compute_spectrogram(.make_sine_psg(), channel = "C3")
  expect_s3_class(sg, "mrpheus_spectrogram")
  expect_equal(nrow(sg$power), 3L)
  expect_equal(ncol(sg$power), length(sg$freqs))
  expect_true(all(sg$freqs >= 0 & sg$freqs <= 40))
})

test_that("compute_spectrogram puts the peak at the input frequency", {
  sg <- compute_spectrogram(.make_sine_psg(freq = 10), channel = "C3", db = FALSE)
  peak <- sg$freqs[which.max(colMeans(sg$power))]
  # nfft = 256 at 100 Hz, so bin width is about 0.39 Hz
  expect_lt(abs(peak - 10), 0.5)
})
