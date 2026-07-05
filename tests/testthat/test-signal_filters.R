# tests/testthat/test-signal_filters.R
#
# Unit tests for signal_filters.R (remove_dc, detect_powerline,
# notch_filter, bandpass_filter) and the preprocess_psg() pipeline.
#
# Approach: synthetic signals with known spectral properties. All expected
# values are derived analytically — no external fixtures needed.

# ── Helpers ──────────────────────────────────────────────────────────────────

# Build a minimal valid mrpheus_psg from scratch.
# Two channels: one EEG (100 Hz) and one EOG (100 Hz).
make_mock_psg <- function(sr = 100, duration_s = 90, seed = 42) {
  set.seed(seed)
  n  <- as.integer(sr * duration_s)
  t  <- seq_len(n) / sr

  # EEG: broadband noise + faint 50 Hz tone for powerline detection
  eeg <- rnorm(n) + 0.5 * sin(2 * pi * 50 * t)
  eog <- rnorm(n) * 0.5

  edf <- list(
    channels = data.frame(
      label       = c("EEG C3", "EOG LOC"),
      sample_rate = c(sr, sr),
      stringsAsFactors = FALSE
    ),
    signals = list(
      `EEG C3`  = list(signal = eeg),
      `EOG LOC` = list(signal = eog)
    )
  )

  cmap <- tibble::tibble(
    label       = c("EEG C3", "EOG LOC"),
    type        = c("EEG",    "EOG"),
    sample_rate = c(sr, sr),
    bad         = c(FALSE, FALSE)
  )

  epoch_s  <- 30L
  ep_samp  <- as.integer(epoch_s * sr)
  n_epochs <- floor(duration_s / epoch_s)

  epochs <- lapply(seq_len(n_epochs), function(i) {
    start <- (i - 1L) * ep_samp + 1L
    list(
      `EEG C3`  = eeg[start:(start + ep_samp - 1L)],
      `EOG LOC` = eog[start:(start + ep_samp - 1L)]
    )
  })

  structure(
    list(
      edf         = edf,
      epochs      = epochs,
      n_epochs    = n_epochs,
      epoch_s     = epoch_s,
      channel_map = cmap
    ),
    class = "mrpheus_psg"
  )
}

# ── remove_dc ────────────────────────────────────────────────────────────────

test_that("remove_dc subtracts the signal mean", {
  x      <- c(10, 12, 11, 13, 9) + 100  # large DC offset
  result <- remove_dc(x)
  expect_equal(mean(result), 0, tolerance = 1e-10)
})

test_that("remove_dc preserves signal length and relative differences", {
  set.seed(1)
  x      <- rnorm(1000) + 50
  result <- remove_dc(x)
  expect_length(result, 1000L)
  expect_equal(diff(result), diff(x), tolerance = 1e-10)
})

test_that("remove_dc is a no-op on a zero-mean signal", {
  x <- c(-1, 1, -1, 1)
  expect_equal(remove_dc(x), x, tolerance = 1e-14)
})

test_that("remove_dc handles NA values without crashing", {
  x <- c(1, 2, NA, 4)
  expect_no_error(remove_dc(x))
})

# ── detect_powerline ─────────────────────────────────────────────────────────

test_that("detect_powerline detects 50 Hz when 50 Hz tone dominates", {
  sr <- 256L
  t  <- seq_len(10L * sr) / sr
  # Strong 50 Hz tone, faint 60 Hz
  sig <- sin(2 * pi * 50 * t) + 0.05 * sin(2 * pi * 60 * t)
  expect_equal(detect_powerline(sig, sr), 50L)
})

test_that("detect_powerline detects 60 Hz when 60 Hz tone dominates", {
  sr <- 256L
  t  <- seq_len(10L * sr) / sr
  # Strong 60 Hz tone, faint 50 Hz
  sig <- sin(2 * pi * 60 * t) + 0.05 * sin(2 * pi * 50 * t)
  expect_equal(detect_powerline(sig, sr), 60L)
})

test_that("detect_powerline returns 50L when sr is too low to see 60 Hz", {
  # Nyquist = 50 Hz, so 60 Hz >= Nyquist → always returns 50
  sr  <- 100L
  sig <- rnorm(10L * sr)
  expect_equal(detect_powerline(sig, sr), 50L)
})

test_that("detect_powerline returns an integer", {
  sr  <- 256L
  sig <- rnorm(10L * sr)
  result <- detect_powerline(sig, sr)
  expect_true(is.integer(result))
  expect_true(result %in% c(50L, 60L))
})

# ── notch_filter ─────────────────────────────────────────────────────────────

test_that("notch_filter preserves signal length", {
  sr  <- 200L
  sig <- rnorm(60L * sr)
  expect_length(notch_filter(sig, sr, freq = 50), length(sig))
})

test_that("notch_filter attenuates the target frequency", {
  sr   <- 256L
  n    <- 30L * sr
  t    <- seq_len(n) / sr
  tone <- sin(2 * pi * 50 * t)
  # SNR: tone + low-power background noise
  sig  <- tone + 0.01 * rnorm(n)

  filt <- notch_filter(sig, sr, freq = 50, harmonics = FALSE)

  # Correlation with the original 50 Hz tone should collapse after notching
  r_before <- abs(cor(sig, tone))
  r_after  <- abs(cor(filt, tone))
  expect_lt(r_after, r_before * 0.1)  # > 90 % of 50 Hz component removed
})

test_that("notch_filter with harmonics removes fundamental and harmonic", {
  sr  <- 256L
  n   <- 30L * sr
  t   <- seq_len(n) / sr
  sig <- sin(2 * pi * 50 * t) + sin(2 * pi * 100 * t) + 0.01 * rnorm(n)

  filt <- notch_filter(sig, sr, freq = 50, harmonics = TRUE)

  tone_100 <- sin(2 * pi * 100 * t)
  r_before <- abs(cor(sig, tone_100))
  r_after  <- abs(cor(filt, tone_100))
  expect_lt(r_after, r_before * 0.1)
})

test_that("notch_filter harmonics = FALSE leaves harmonics intact", {
  sr  <- 256L
  n   <- 30L * sr
  t   <- seq_len(n) / sr
  tone_100 <- sin(2 * pi * 100 * t)
  sig      <- sin(2 * pi * 50 * t) + tone_100

  filt <- notch_filter(sig, sr, freq = 50, harmonics = FALSE)

  # 100 Hz component should survive
  expect_gt(abs(cor(filt, tone_100)), 0.9)
})

test_that("notch_filter returns a numeric vector", {
  sig  <- rnorm(1000)
  filt <- notch_filter(sig, sr = 200, freq = 50, harmonics = FALSE)
  expect_true(is.numeric(filt))
})

# ── bandpass_filter ───────────────────────────────────────────────────────────

test_that("bandpass_filter preserves signal length", {
  sr  <- 200L
  sig <- rnorm(60L * sr)
  expect_length(bandpass_filter(sig, sr, 0.3, 35), length(sig))
})

test_that("bandpass_filter preserves in-band signal and attenuates out-of-band", {
  sr <- 256L
  n  <- 20L * sr
  t  <- seq_len(n) / sr

  inband   <- sin(2 * pi * 10 * t)   # 10 Hz — in [0.3, 35] Hz passband
  outband  <- sin(2 * pi * 80 * t)   # 80 Hz — above passband
  mixed    <- inband + outband

  filt <- bandpass_filter(mixed, sr, low_hz = 0.3, high_hz = 35)

  # In-band: high correlation with original 10 Hz component.
  # Threshold is 0.95 (not 1.0) to allow for IIR filter edge effects at the
  # start/end of a 20 s signal.
  expect_gt(cor(filt, inband), 0.95)

  # Out-of-band attenuation: correlation of the filtered signal with the
  # original 80 Hz tone should be well below the in-band correlation.
  # Using correlation (not RMS residual) avoids conflating passband ripple
  # with the out-of-band component.
  r_inband  <- cor(filt, inband)
  r_outband <- abs(cor(filt, outband))
  expect_lt(r_outband, r_inband * 0.5)
})

test_that("bandpass_filter with low_hz = 0 acts as a low-pass filter", {
  sr  <- 200L
  n   <- 20L * sr
  t   <- seq_len(n) / sr
  sig <- sin(2 * pi * 5 * t) + sin(2 * pi * 60 * t)

  filt     <- bandpass_filter(sig, sr, low_hz = 0, high_hz = 30)
  tone_60  <- sin(2 * pi * 60 * t)
  expect_lt(sqrt(mean((filt - sin(2 * pi * 5 * t))^2)), 0.1)
})

test_that("bandpass_filter clamps high_hz at Nyquist and acts as high-pass", {
  sr   <- 200L
  n    <- 20L * sr
  t    <- seq_len(n) / sr
  sig  <- sin(2 * pi * 0.05 * t) + sin(2 * pi * 10 * t)

  # high_hz > Nyquist → should fall back to high-pass at low_hz
  filt    <- bandpass_filter(sig, sr, low_hz = 1, high_hz = 200)
  tone_lf <- sin(2 * pi * 0.05 * t)
  # Low-frequency component should be heavily attenuated
  expect_lt(abs(cor(filt, tone_lf)), 0.3)
})

test_that("bandpass_filter returns a numeric vector", {
  sig  <- rnorm(1000)
  filt <- bandpass_filter(sig, sr = 200, low_hz = 0.3, high_hz = 35)
  expect_true(is.numeric(filt))
})

# ── preprocess_psg ────────────────────────────────────────────────────────────

test_that("preprocess_psg returns an mrpheus_psg object", {
  psg  <- make_mock_psg()
  out  <- preprocess_psg(psg, powerline_freq = 50L, verbose = FALSE)
  expect_s3_class(out, "mrpheus_psg")
})

test_that("preprocess_psg preserves epoch count and duration", {
  psg <- make_mock_psg()
  out <- preprocess_psg(psg, powerline_freq = 50L, verbose = FALSE)
  expect_equal(out$n_epochs, psg$n_epochs)
  expect_equal(out$epoch_s,  psg$epoch_s)
})

test_that("preprocess_psg preserves epoch signal lengths", {
  psg <- make_mock_psg()
  out <- preprocess_psg(psg, powerline_freq = 50L, verbose = FALSE)
  # Every epoch and channel should have the same sample count as the input
  for (i in seq_len(out$n_epochs)) {
    for (ch in out$channel_map$label) {
      expect_equal(length(out$epochs[[i]][[ch]]),
                   length(psg$epochs[[i]][[ch]]))
    }
  }
})

test_that("preprocess_psg channel_rename updates channel_map and signals", {
  psg <- make_mock_psg()
  out <- preprocess_psg(
    psg,
    channel_rename = c("EEG C3" = "C3"),
    powerline_freq = 50L,
    verbose        = FALSE
  )
  expect_true("C3" %in% out$channel_map$label)
  expect_false("EEG C3" %in% out$channel_map$label)
  expect_true("C3" %in% names(out$edf$signals))
})

test_that("preprocess_psg channel_rename emits a message on unknown labels", {
  # cli::cli_alert_warning emits via message(), not warning(), so expect_message.
  psg <- make_mock_psg()
  expect_message(
    preprocess_psg(
      psg,
      channel_rename = c("NONEXISTENT" = "X"),
      powerline_freq = 50L,
      verbose        = TRUE
    ),
    regexp = "NONEXISTENT"
  )
})

test_that("preprocess_psg with dc = FALSE skips DC removal", {
  # A signal with a known offset; after dc=FALSE it should retain offset in
  # continuous signal (DC is preserved before bandpass removes it anyway).
  # We simply verify it runs without error.
  psg <- make_mock_psg()
  expect_no_error(
    preprocess_psg(psg, dc = FALSE, powerline_freq = 50L, verbose = FALSE)
  )
})

test_that("preprocess_psg auto-detects powerline when powerline_freq = NULL", {
  psg <- make_mock_psg()   # mock has faint 50 Hz tone in EEG channel
  expect_no_error(
    preprocess_psg(psg, powerline_freq = NULL, verbose = FALSE)
  )
})

test_that("preprocess_psg stops on non-mrpheus_psg input", {
  expect_error(preprocess_psg(list()), regexp = "mrpheus_psg")
})
