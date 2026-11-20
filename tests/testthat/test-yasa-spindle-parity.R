# tests/testthat/test-yasa-spindle-parity.R
#
# compute_spindles() with its default arguments must reproduce
# yasa.spindles_detect() (YASA 0.8.1, default arguments): the same events, at the
# same sample positions, with the same properties. Fixtures
# (data-raw/generate_yasa_spindle_parity_fixture.py): the events YASA finds in a
# deterministic synthetic EEG with spindle bursts at 100 Hz and at 256 Hz (12
# epochs, three of them excluded by the stage mask) and two variants at 100 Hz
# ("loud": strong sigma activity in an excluded epoch; "x4": the signal times 4),
# and, at 100 Hz, YASA's intermediate features (relative sigma power of the STFT,
# moving correlation and moving RMS). The signal is generated here from the same
# formulas, so no signal files are stored.
#
# Each of these was put back to what it was before (or removed), one at a time,
# and the tests fail:
#   the broadband and the sigma FIR filters (different transition bands)
#   the STFT relative power, the moving correlation and the moving RMS criteria
#   the RMS threshold (trimmed standard deviation, mean within the stage mask,
#     cap at 10 uV)
#   the 0.1 s smoothing of the criteria and its alignment
#   the merging of events closer than 500 ms
#   the duration limits (exclusive), and the stage mask
#   the properties of an event (detrending, relative power, oscillation count)
# Agreement with YASA is exact in time (same samples) and about 1e-11 in the
# properties.

.sp_lcg <- function(n, seed) {       # Park-Miller; 16807 * state < 2^53, exact in doubles
  out <- numeric(n); state <- seed
  for (i in seq_len(n)) {
    state <- (16807 * state) %% 2147483647
    out[i] <- state / 2147483647 - 0.5
  }
  out
}

# (centre s, half-width s, amplitude uV, frequency Hz)
.sp_bursts <- matrix(c(
  12.0, 0.45, 38, 13.5,   31.0, 0.60, 34, 12.8,   47.5, 0.35, 42, 14.2,   48.4, 0.30, 36, 13.0,
  75.0, 0.55, 40, 13.6,  101.0, 0.70, 33, 13.2,  102.6, 0.40, 35, 14.0,  128.0, 0.50, 45, 13.8,
  151.0, 0.25, 40, 13.5, 166.0, 0.80, 37, 12.6,  196.0, 0.55, 41, 13.9,  222.0, 0.45, 39, 13.3,
  241.0, 0.65, 36, 14.4, 263.0, 0.50, 44, 13.1,  283.0, 0.60, 38, 13.7,  305.0, 0.40, 43, 13.4,
  329.0, 0.55, 35, 12.9, 347.0, 0.70, 40, 13.6), ncol = 4, byrow = TRUE)

.sp_signal <- function(fs, n_epochs = 12L) {
  n  <- fs * 30L * n_epochs
  t  <- (seq_len(n) - 1L) / fs
  tp <- 2 * pi
  x  <- 28 * (1 + 0.4 * sin(tp * t / 150)) * sin(tp * 0.8 * t + 0.4) +
        10 * sin(tp * 6.1 * t + 1.0) + 14 * .sp_lcg(n, 12345)
  for (k in seq_len(nrow(.sp_bursts))) {
    b <- .sp_bursts[k, ]
    x <- x + b[3] * exp(-((t - b[1]) / b[2])^2) * sin(tp * b[4] * (t - b[1]))
  }
  x
}

.sp_psg <- function(fs, n_epochs = 12L) {
  sig <- .sp_signal(fs, n_epochs)
  spe <- fs * 30L
  structure(
    list(
      epochs = lapply(seq_len(n_epochs), function(i)
        list(EEG = sig[((i - 1L) * spe + 1L):(i * spe)])),
      n_epochs = n_epochs, epoch_s = 30,
      channel_map = tibble::tibble(label = "EEG", type = "EEG", sample_rate = fs, bad = FALSE)
    ),
    class = "mrpheus_psg"
  )
}

.sp_stages <- c(1:3, 5:7, 9:11)       # epochs 4, 8 and 12 are excluded

test_that("the synthetic EEG matches the one the spindle fixture was made from", {
  # printed by data-raw/generate_yasa_spindle_parity_fixture.py
  expect_equal(.sp_signal(100)[1:2], c(13.67105483182771, 26.68838660591996), tolerance = 1e-10)
  expect_equal(.sp_signal(256)[1:2], c(13.67105483182771, 25.21111149056969), tolerance = 1e-10)
})

test_that("the STFT relative power and the moving correlation and RMS equal YASA's", {
  fix <- utils::read.csv(testthat::test_path("fixtures", "yasa_sp_features.csv"))
  fs  <- 100
  x   <- .sp_signal(fs)
  tb  <- mrpheus:::.mne_auto_trans(fs, 1, 30)
  xb  <- mrpheus:::.mne_fir_apply(x, mrpheus:::.mne_fir_bandpass_design(fs, 1, 30, tb[["l"]], tb[["h"]]))
  xs  <- mrpheus:::.mne_fir_apply(x, mrpheus:::.mne_fir_bandpass_design(fs, 12, 15, 1.5, 1.5))

  st <- mrpheus:::.stft_rel_pow(xb, fs, c(12, 15), c(1, 30))
  r  <- fix[fix$kind == "rel_pow", ]
  expect_equal(st$t[r$index],   r$time,  tolerance = 1e-12)
  expect_equal(st$rel[r$index], r$value, tolerance = 1e-9)

  mc <- mrpheus:::.moving_transform(xs, xb, fs, "corr")
  r  <- fix[fix$kind == "corr", ]
  expect_equal(mc[r$index], r$value, tolerance = 1e-9)

  mr <- mrpheus:::.moving_transform(xs, NULL, fs, "rms")
  r  <- fix[fix$kind == "rms", ]
  expect_equal(mr[r$index], r$value, tolerance = 1e-9)
})

# the signal of a fixture variant (see the generator)
.sp_variant_psg <- function(variant, fs) {
  psg <- .sp_psg(fs)
  spe <- fs * 30L
  if (variant == "loud") {
    t   <- (seq_len(spe) - 1L) / fs + 90                       # epoch 4: 90-120 s
    psg$epochs[[4]]$EEG <- psg$epochs[[4]]$EEG +
      120 * sin(2 * pi * 13.5 * t) * (1 + 0.5 * sin(2 * pi * 0.3 * t))
  } else if (variant == "x4") {
    psg$epochs <- lapply(psg$epochs, function(e) list(EEG = 4 * e$EEG))
  }
  psg
}

test_that("compute_spindles finds the events of yasa.spindles_detect", {
  fix <- utils::read.csv(testthat::test_path("fixtures", "yasa_sp_events.csv"))

  for (key in unique(paste(fix$variant, fix$fs))) {
    y   <- fix[paste(fix$variant, fix$fs) == key, ]
    psg <- .sp_variant_psg(y$variant[1], y$fs[1])
    sp  <- suppressMessages(compute_spindles(psg, channel = "EEG", stages = .sp_stages))
    info <- key

    expect_equal(nrow(sp), nrow(y), info = info)
    # YASA's times are seconds from the first sample; mrpheus gives start and end
    # relative to the epoch of the first sample
    expect_equal((sp$epoch - 1) * 30 + sp$start_s, y$Start, tolerance = 1e-9, info = paste(info, "start"))
    expect_equal((sp$epoch - 1) * 30 + sp$end_s,   y$End,   tolerance = 1e-9, info = paste(info, "end"))
    expect_equal(sp$duration_s,     y$Duration,     tolerance = 1e-9, info = paste(info, "duration"))
    expect_equal(sp$peak_s,         y$Peak,         tolerance = 1e-9, info = paste(info, "peak"))
    expect_equal(sp$amplitude_uv,   y$Amplitude,    tolerance = 1e-8, info = paste(info, "amplitude"))
    expect_equal(sp$amp_filtered_uv, y$AmpFiltered, tolerance = 1e-8, info = paste(info, "filtered amplitude"))
    expect_equal(sp$rms_uv,         y$RMS,          tolerance = 1e-8, info = paste(info, "rms"))
    expect_equal(sp$abs_power,      y$AbsPower,     tolerance = 1e-8, info = paste(info, "abs power"))
    expect_equal(sp$rel_power,      y$RelPower,     tolerance = 1e-8, info = paste(info, "rel power"))
    expect_equal(sp$peak_freq_hz,   y$Frequency,    tolerance = 1e-8, info = paste(info, "frequency"))
    expect_equal(sp$oscillations,   as.integer(y$Oscillations), info = paste(info, "oscillations"))
    expect_equal(sp$symmetry,       y$Symmetry,     tolerance = 1e-8, info = paste(info, "symmetry"))
  }
})

test_that("the duration limits are exclusive", {
  psg <- .sp_psg(100)
  d0  <- suppressMessages(compute_spindles(psg, stages = .sp_stages))$duration_s[1]
  lo  <- suppressMessages(compute_spindles(psg, stages = .sp_stages, duration_s = c(d0, 2)))
  hi  <- suppressMessages(compute_spindles(psg, stages = .sp_stages, duration_s = c(0.5, d0)))
  expect_false(any(abs(lo$duration_s - d0) < 1e-9))
  expect_false(any(abs(hi$duration_s - d0) < 1e-9))
})

test_that("the fixture exercises merging, the stage mask and the duration limits", {
  fix <- utils::read.csv(testthat::test_path("fixtures", "yasa_sp_events.csv"))
  fix <- fix[fix$variant == "base", ]
  for (fs in unique(fix$fs)) {
    y <- fix[fix$fs == fs, ]
    expect_gte(nrow(y), 10L)
    expect_gt(max(y$Duration), 1.5)                  # the pair of bursts 0.9 s apart is merged
    # nothing starts in the excluded epochs 4, 8 and 12 (90-120, 210-240, 330-360 s)
    expect_false(any((floor(y$Start / 30) + 1) %in% c(4, 8, 12)))
    all <- suppressMessages(compute_spindles(.sp_psg(fs), channel = "EEG"))
    expect_gt(nrow(all), nrow(y))
  }
})
