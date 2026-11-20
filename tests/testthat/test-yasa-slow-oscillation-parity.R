# tests/testthat/test-yasa-slow-oscillation-parity.R
#
# compute_slow_oscillations() with its default arguments must reproduce
# yasa.sw_detect() (YASA 0.8.1, default arguments): the same events, at the same
# sample positions, with the same amplitudes. Fixtures
# (data-raw/generate_yasa_slow_oscillation_parity_fixture.py): the events YASA
# finds in a deterministic synthetic EEG at 100 Hz and at 256 Hz (10 epochs,
# two of them excluded by the stage mask, so that events are cut by the mask and
# by epoch boundaries), and taps of the FIR filter that sw_detect applies. The
# synthetic EEG is generated here from the same formulas. The crafted signal of
# the last test (and the events YASA finds in it) is stored as it is.
#
# Each of these was put back to what it was before (or removed), one at a time,
# and the tests fail:
#   YASA's removal of events that share a start with another event
#   the FIR band-pass (a different transition band or window, or a per-epoch filter)
#   pairing each negative peak with the closest following positive peak
#   the trough and peak amplitude limits
#   the stage mask applying to the peaks
#   limits that are exclusive: peak-to-peak amplitude, duration of each phase
#   the check that the duration equals the sum of the two phases
#   zero crossings defined by x > 0 (an exact zero counts as not positive)
# Two deviations cannot be detected because they change nothing: removing only
# the shared-end check (when the duration rule holds, a shared end implies a
# shared start), and the wrap-around of the previous zero crossing before the
# first one (the event is rejected either way).
# Agreement with YASA is exact in time (same samples) and about 1e-12 in amplitude.

.sw_lcg <- function(n, seed) {       # Park-Miller; 16807 * state < 2^53, exact in doubles
  out <- numeric(n); state <- seed
  for (i in seq_len(n)) {
    state <- (16807 * state) %% 2147483647
    out[i] <- state / 2147483647 - 0.5
  }
  out
}

.sw_signal <- function(fs, n_epochs = 10L) {
  n  <- fs * 30L * n_epochs
  t  <- (seq_len(n) - 1L) / fs
  tp <- 2 * pi
  env <- 1 + 0.5 * sin(tp * t / 180)
  70 * env * (sin(tp * 0.75 * t) + 0.45 * sin(tp * 1.15 * t + 0.5)) +
    14 * sin(tp * 0.3 * t + 2) + 18 * .sw_lcg(n, 12345)
}

.sw_psg <- function(fs, n_epochs = 10L) {
  sig <- .sw_signal(fs, n_epochs)
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

.sw_stages <- c(1:3, 5:7, 9:10)       # epochs 4 and 8 are excluded

test_that("the synthetic EEG matches the one the slow oscillation fixture was made from", {
  # printed by data-raw/generate_yasa_slow_oscillation_parity_fixture.py
  expect_equal(.sw_signal(100)[1:2], c(20.5711659547289, 38.989209797359734),  tolerance = 1e-10)
  expect_equal(.sw_signal(256)[1:2], c(20.5711659547289, 35.864506206926606), tolerance = 1e-10)
})

test_that("the FIR band-pass design equals mne.filter.create_filter", {
  fix <- utils::read.csv(testthat::test_path("fixtures", "yasa_sw_filter.csv"))
  for (fs in unique(fix$fs)) {
    f <- fix[fix$fs == fs, ]
    h <- mrpheus:::.mne_fir_bandpass_design(fs, 0.3, 1.5, 0.2, 0.2)
    expect_length(h, f$n_taps[1])
    expect_equal(h[f$tap], f$h, tolerance = 1e-12, info = paste("fs =", fs))
    expect_equal(h, rev(h))                          # symmetric: linear phase
  }
})

test_that("compute_slow_oscillations finds the events of yasa.sw_detect", {
  fix <- utils::read.csv(testthat::test_path("fixtures", "yasa_sw_events.csv"))

  for (fs in unique(fix$fs)) {
    y   <- fix[fix$fs == fs, ]
    psg <- .sw_psg(fs)
    so  <- suppressMessages(compute_slow_oscillations(psg, channel = "EEG", stages = .sw_stages))
    info <- paste("fs =", fs)

    expect_equal(nrow(so), nrow(y), info = info)
    # YASA's times are seconds from the first sample; mrpheus gives start and end
    # relative to the epoch of the negative peak
    start <- (so$epoch - 1) * 30 + so$start_s
    end   <- (so$epoch - 1) * 30 + so$end_s
    expect_equal(start,              y$Start,       tolerance = 1e-9, info = paste(info, "start"))
    expect_equal(so$neg_peak_s,      y$NegPeak,     tolerance = 1e-9, info = paste(info, "neg peak"))
    expect_equal(so$mid_crossing_s,  y$MidCrossing, tolerance = 1e-9, info = paste(info, "mid crossing"))
    expect_equal(so$pos_peak_s,      y$PosPeak,     tolerance = 1e-9, info = paste(info, "pos peak"))
    expect_equal(end,                y$End,         tolerance = 1e-9, info = paste(info, "end"))
    expect_equal(so$duration_s,      y$Duration,    tolerance = 1e-9, info = paste(info, "duration"))
    expect_equal(so$neg_peak_uv,     y$ValNegPeak,  tolerance = 1e-8, info = paste(info, "neg value"))
    expect_equal(so$pos_peak_uv,     y$ValPosPeak,  tolerance = 1e-8, info = paste(info, "pos value"))
    expect_equal(so$ptp_uv,          y$PTP,         tolerance = 1e-8, info = paste(info, "ptp"))
    expect_equal(so$slope_uv_s,      y$Slope,       tolerance = 1e-8, info = paste(info, "slope"))
  }
})

test_that("the fixture exercises the stage mask and events that cross epoch boundaries", {
  fix <- utils::read.csv(testthat::test_path("fixtures", "yasa_sw_events.csv"))
  for (fs in unique(fix$fs)) {
    y <- fix[fix$fs == fs, ]
    expect_gt(nrow(y), 100L)
    crossing <- floor(y$Start / 30) != floor(y$NegPeak / 30) | floor(y$End / 30) != floor(y$NegPeak / 30)
    expect_gte(sum(crossing), 3L)
    # no negative peak in an excluded epoch, and the mask removes events
    expect_false(any((floor(y$NegPeak / 30) + 1) %in% c(4, 8)))
    all <- suppressMessages(compute_slow_oscillations(.sw_psg(fs), channel = "EEG"))
    expect_gt(nrow(all), nrow(y))
  }
})

test_that("the detection equals YASA's on a crafted signal that reaches its edge cases", {
  # An already band-passed signal (integer uV, 100 Hz) with peak plateaus,
  # double-humped waves, durations and amplitudes exactly on a limit, exact
  # zeros between waves, a wave before the first zero crossing and a stage mask,
  # and the events of YASA's own detection code (sw_detect with its band-pass
  # replaced by the identity).
  sig <- utils::read.csv(testthat::test_path("fixtures", "yasa_sw_crafted_signal.csv"))
  y   <- utils::read.csv(testthat::test_path("fixtures", "yasa_sw_crafted_events.csv"))
  ev  <- mrpheus:::.sw_events(sig$x, 100, sig$mask == 1,
                              dur_neg = c(0.3, 1.5), dur_pos = c(0.1, 1.0),
                              amp_neg = c(40, 200), amp_pos = c(10, 150),
                              amp_ptp = c(75, 350))
  expect_equal(nrow(ev), nrow(y))
  for (col in c("Start", "NegPeak", "MidCrossing", "PosPeak", "End", "Duration"))
    expect_equal(ev[[col]], y[[col]], tolerance = 1e-9, info = col)
  for (col in c("ValNegPeak", "ValPosPeak", "PTP", "Slope", "Frequency"))
    expect_equal(ev[[col]], y[[col]], tolerance = 1e-9, info = col)
})
