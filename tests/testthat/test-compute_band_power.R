# tests/testthat/test-compute_band_power.R
#
# Unit tests for compute_band_power() and the internal .welch_psd_bp().
#
# Approach: synthetic signals with known spectral properties so expected
# values can be derived analytically without external fixtures.

# ── Helpers ──────────────────────────────────────────────────────────────────

make_bp_psg <- function(sr = 100, epoch_s = 30, n_epochs = 3, seed = 1) {
  set.seed(seed)
  n_ep  <- as.integer(sr * epoch_s)
  n_tot <- n_ep * n_epochs
  sig   <- rnorm(n_tot)

  edf <- list(
    channels = data.frame(label = "EEG C3", sample_rate = sr,
                          stringsAsFactors = FALSE),
    signals  = list(`EEG C3` = list(signal = sig))
  )
  cmap <- tibble::tibble(
    label = "EEG C3", type = "EEG", sample_rate = sr, bad = FALSE
  )
  epochs <- lapply(seq_len(n_epochs), function(i) {
    start <- (i - 1L) * n_ep + 1L
    list(`EEG C3` = sig[start:(start + n_ep - 1L)])
  })
  structure(
    list(edf = edf, epochs = epochs, n_epochs = n_epochs,
         epoch_s = epoch_s, channel_map = cmap),
    class = "mrpheus_psg"
  )
}

# ── .welch_psd_bp ─────────────────────────────────────────────────────────────

test_that(".welch_psd_bp returns correct frequency vector length", {
  nfft  <- 1024L
  n_one <- nfft %/% 2L + 1L
  set.seed(1)
  sig   <- rnorm(3000)
  psd   <- mrpheus:::.welch_psd_bp(sig, fs = 100, win_sec = 4,
                                    noverlap = 200L, nfft = nfft)
  expect_length(psd$freq, n_one)
  expect_length(psd$spec, n_one)
})

test_that(".welch_psd_bp frequency vector starts at 0 and ends at Nyquist", {
  set.seed(1)
  psd <- mrpheus:::.welch_psd_bp(rnorm(3000), fs = 100,
                                  win_sec = 4, noverlap = 200L, nfft = 1024L)
  expect_equal(psd$freq[1], 0)
  expect_equal(psd$freq[length(psd$freq)], 50, tolerance = 1e-10)
})

test_that(".welch_psd_bp PSD is non-negative", {
  set.seed(2)
  psd <- mrpheus:::.welch_psd_bp(rnorm(3000), fs = 100,
                                  win_sec = 4, noverlap = 200L, nfft = 1024L)
  expect_true(all(psd$spec >= 0))
})

test_that(".welch_psd_bp detects tonal power at known frequency", {
  sr  <- 100L
  t   <- seq_len(30L * sr) / sr
  sig <- sin(2 * pi * 10 * t)    # pure 10 Hz tone

  psd      <- mrpheus:::.welch_psd_bp(sig, fs = sr, win_sec = 4,
                                       noverlap = 200L, nfft = 1024L)
  peak_idx <- which.max(psd$spec)
  expect_equal(psd$freq[peak_idx], 10, tolerance = 0.2)  # peak near 10 Hz
})

test_that(".welch_psd_bp raises nfft silently when nfft < nperseg", {
  set.seed(3)
  sig <- rnorm(3000)
  # nfft=128 < nperseg=400; should not error
  expect_no_error(
    mrpheus:::.welch_psd_bp(sig, fs = 100, win_sec = 4,
                             noverlap = 50L, nfft = 128L)
  )
})

# ── compute_band_power ────────────────────────────────────────────────────────

test_that("compute_band_power returns a tibble with correct columns", {
  psg <- make_bp_psg()
  bp  <- compute_band_power(psg)

  expect_s3_class(bp, "tbl_df")
  expect_true(all(c("epoch", "channel", "delta", "theta", "alpha",
                    "sigma", "beta", "gamma", "total_power") %in% names(bp)))
})

test_that("compute_band_power returns one row per epoch per channel", {
  psg <- make_bp_psg(n_epochs = 3)
  bp  <- compute_band_power(psg)
  expect_equal(nrow(bp), 3L)   # 3 epochs x 1 channel
})

test_that("compute_band_power epoch column runs 1 to n_epochs", {
  psg <- make_bp_psg(n_epochs = 4)
  bp  <- compute_band_power(psg)
  expect_equal(sort(unique(bp$epoch)), 1:4)
})

test_that("compute_band_power band powers are non-negative", {
  psg       <- make_bp_psg()
  bp        <- compute_band_power(psg)
  band_cols <- c("delta", "theta", "alpha", "sigma", "beta", "gamma")
  expect_true(all(bp[band_cols] >= 0, na.rm = TRUE))
})

test_that("compute_band_power relative = TRUE sums bands to 1", {
  psg <- make_bp_psg()
  bp  <- compute_band_power(psg, relative = TRUE)
  band_cols <- c("delta", "theta", "alpha", "sigma", "beta", "gamma")
  band_sums <- rowSums(bp[band_cols], na.rm = TRUE)
  expect_equal(band_sums, rep(1, nrow(bp)), tolerance = 1e-10)
})

test_that("compute_band_power relative = FALSE preserves absolute units", {
  psg <- make_bp_psg()
  bp  <- compute_band_power(psg, relative = FALSE)
  # total_power should equal sum of bands
  band_cols <- c("delta", "theta", "alpha", "sigma", "beta", "gamma")
  band_sums <- rowSums(bp[band_cols], na.rm = TRUE)
  expect_equal(bp$total_power, band_sums, tolerance = 1e-10)
})

test_that("compute_band_power custom bands are respected", {
  psg <- make_bp_psg()
  bp  <- compute_band_power(
    psg,
    bands   = list(slow = c(0.5, 2), fast = c(20, 40))
  )
  expect_named(bp[-(1:2)], c("slow", "fast", "total_power"))
})

test_that("compute_band_power channels argument subsets correctly", {
  # Two-channel PSG
  sr <- 100L; n <- 3000L
  set.seed(5)
  edf <- list(
    channels = data.frame(label = c("EEG C3", "EEG C4"), sample_rate = c(sr, sr),
                          stringsAsFactors = FALSE),
    signals  = list(`EEG C3` = list(signal = rnorm(n)),
                    `EEG C4` = list(signal = rnorm(n)))
  )
  cmap <- tibble::tibble(
    label = c("EEG C3", "EEG C4"), type = "EEG",
    sample_rate = sr, bad = FALSE
  )
  ep_n <- as.integer(30 * sr)
  epochs <- lapply(1L, function(i) list(
    `EEG C3` = rnorm(ep_n), `EEG C4` = rnorm(ep_n)
  ))
  psg2 <- structure(
    list(edf = edf, epochs = epochs, n_epochs = 1L,
         epoch_s = 30, channel_map = cmap),
    class = "mrpheus_psg"
  )

  bp_all  <- compute_band_power(psg2)
  bp_one  <- compute_band_power(psg2, channels = "EEG C3")

  expect_equal(nrow(bp_all), 2L)   # 1 epoch x 2 channels
  expect_equal(nrow(bp_one), 1L)
  expect_equal(bp_one$channel, "EEG C3")
})

test_that("compute_band_power stops on non-mrpheus_psg input", {
  expect_error(compute_band_power(list()), regexp = "mrpheus_psg")
})

test_that("compute_band_power delta power is higher for a delta-band signal", {
  sr    <- 100L
  n_ep  <- as.integer(30 * sr)
  t     <- seq_len(n_ep) / sr
  # Signal dominated by 2 Hz (delta band)
  delta_sig <- sin(2 * pi * 2 * t)

  edf <- list(
    channels = data.frame(label = "EEG C3", sample_rate = sr,
                          stringsAsFactors = FALSE),
    signals  = list(`EEG C3` = list(signal = delta_sig))
  )
  cmap <- tibble::tibble(label = "EEG C3", type = "EEG",
                          sample_rate = sr, bad = FALSE)
  psg_delta <- structure(
    list(edf = edf, epochs = list(list(`EEG C3` = delta_sig)),
         n_epochs = 1L, epoch_s = 30, channel_map = cmap),
    class = "mrpheus_psg"
  )

  bp <- compute_band_power(psg_delta, relative = TRUE)
  # Delta should be the dominant band
  band_cols <- c("delta", "theta", "alpha", "sigma", "beta", "gamma")
  expect_equal(which.max(unlist(bp[band_cols])), 1L)  # delta is column 1
})
