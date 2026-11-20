# tests/testthat/test-compute_spindles.R
#
# Behaviour of compute_spindles() and its helpers. Agreement with
# yasa.spindles_detect() is tested in test-yasa-spindle-parity.R. Expected values
# of the helpers come from scipy and YASA 0.8.1.
#
# All fixtures are fully synthetic (no participant-derived data).

.make_sp_psg <- function(n_epochs = 3L, sr = 100, burst_s = c(10, 40, 70), epoch_s = 30, amp = 40) {
  n  <- as.integer(epoch_s * sr * n_epochs)
  tt <- (seq_len(n) - 1L) / sr
  set.seed(3)
  x  <- 8 * rnorm(n)
  for (c in burst_s) x <- x + amp * exp(-((tt - c) / 0.5)^2) * sin(2 * pi * 13.5 * (tt - c))
  spe <- as.integer(epoch_s * sr)
  structure(
    list(
      epochs      = lapply(seq_len(n_epochs), function(i) list(C3 = x[((i - 1L) * spe + 1L):(i * spe)])),
      n_epochs    = n_epochs,
      epoch_s     = epoch_s,
      channel_map = tibble::tibble(label = "C3", type = "EEG", sample_rate = sr, bad = FALSE)
    ),
    class = "mrpheus_psg"
  )
}

test_that("compute_spindles finds planted spindles and returns the documented columns", {
  sp <- compute_spindles(.make_sp_psg(), channel = "C3")
  expect_s3_class(sp, "tbl_df")
  expect_gte(nrow(sp), 2L)
  expect_named(sp, c("epoch", "start_s", "end_s", "duration_s", "peak_freq_hz", "rms_uv", "channel",
                     "peak_s", "amplitude_uv", "amp_filtered_uv", "abs_power", "rel_power",
                     "oscillations", "symmetry"))
  expect_true(all(sp$duration_s > 0.5 & sp$duration_s < 2))
  expect_true(all(sp$peak_freq_hz > 12 & sp$peak_freq_hz < 15))
  expect_true(all(sp$symmetry >= 0 & sp$symmetry <= 1))
  expect_true(all(sp$rel_power > 0 & sp$rel_power <= 1))
  expect_true(all(sp$epoch %in% 1:3))
})

test_that("stages restricts the epochs in which spindles can occur", {
  psg <- .make_sp_psg()
  all <- suppressMessages(compute_spindles(psg))
  one <- suppressMessages(compute_spindles(psg, stages = 2L))
  expect_gt(nrow(all), nrow(one))
  expect_true(all(one$epoch == 2L))
})

test_that("a threshold set to NULL switches that criterion off, and at least one is needed", {
  psg <- .make_sp_psg()
  base <- suppressMessages(compute_spindles(psg))
  expect_no_error(suppressMessages(compute_spindles(psg, thresh_corr = NULL)))
  expect_no_error(suppressMessages(compute_spindles(psg, thresh_rel_pow = NULL, thresh_rms = NULL)))
  expect_error(compute_spindles(psg, thresh_rel_pow = NULL, thresh_corr = NULL, thresh_rms = NULL),
               "At least one")
  # a stricter criterion cannot add events
  expect_lte(nrow(suppressMessages(compute_spindles(psg, thresh_rel_pow = 0.9))), nrow(base))
})

test_that("a channel in the wrong unit is skipped with a warning", {
  psg <- .make_sp_psg()
  psg$epochs <- lapply(psg$epochs, function(e) list(C3 = e$C3 * 1e-6))   # volts, not uV
  expect_warning(sp <- compute_spindles(psg), "Wrong data amplitude")
  expect_equal(nrow(sp), 0L)
  expect_true("symmetry" %in% names(sp))
})

test_that("a channel without sigma activity returns an empty result", {
  psg <- .make_sp_psg()
  tt  <- (seq_len(3000L) - 1L) / 100
  psg$epochs <- lapply(psg$epochs, function(e) list(C3 = 30 * sin(2 * pi * 6 * tt)))   # theta only
  sp <- suppressMessages(compute_spindles(psg))
  expect_equal(nrow(sp), 0L)
  expect_named(sp, c("epoch", "start_s", "end_s", "duration_s", "peak_freq_hz", "rms_uv", "channel",
                     "peak_s", "amplitude_uv", "amp_filtered_uv", "abs_power", "rel_power",
                     "oscillations", "symmetry"))
})

test_that(".sp_rms_threshold is the mean plus trimmed SD within the mask, at most 10 uV", {
  # expected values from YASA's formula: x[mask].mean() + 1.5 * trimbothstd(x[mask], cut = 0.10), min(., 10)
  x <- seq(1, 3, length.out = 100)
  m <- c(rep(TRUE, 80), rep(FALSE, 20))
  expect_equal(mrpheus:::.sp_rms_threshold(x, m, 1.5), 2.3621915169199577, tolerance = 1e-12)
  # values outside the mask do not matter
  z <- c(seq(1, 3, length.out = 80), rep(500, 20))
  expect_equal(mrpheus:::.sp_rms_threshold(z, m, 1.5), 2.7070501287984277, tolerance = 1e-12)
  # 85.357 without the cap
  expect_equal(mrpheus:::.sp_rms_threshold(as.numeric(1:100), rep(TRUE, 100), 1.5), 10)
})

test_that(".merge_close fills gaps shorter than the minimum distance (YASA _merge_close)", {
  # gaps of 8 and 19 samples are filled at 100 Hz and 500 ms (50 samples); the gap of 159 is not
  expect_equal(mrpheus:::.merge_close(c(10L, 11L, 12L, 20L, 21L, 40L, 41L, 200L), 500, 100),
               c(10:41, 200L))
  expect_equal(mrpheus:::.merge_close(c(1L, 2L, 3L), 500, 100), c(1L, 2L, 3L))
})

test_that(".detrend_linear equals scipy.signal.detrend(type = 'linear')", {
  expect_equal(mrpheus:::.detrend_linear(c(1, 3, 2, 5, 4, 8)),
               c(0.095238095238, 0.92380952381, -1.247619047619, 0.580952380952,
                 -1.590476190476, 1.238095238095), tolerance = 1e-10)
})

test_that(".peaks_prominence follows scipy find_peaks(distance, prominence)", {
  x <- c(0, 3, 1, 2, 0.5, 4, 0, 1, 0.2, 2.5, 0.1)
  a <- mrpheus:::.peaks_prominence(x, 3)
  expect_equal(a$peaks, c(2L, 6L, 10L))
  expect_equal(a$prom, c(2.5, 4, 2.4))
  b <- mrpheus:::.peaks_prominence(x, 2)
  expect_equal(b$peaks, c(2L, 4L, 6L, 8L, 10L))
  expect_equal(b$prom, c(2.5, 1, 4, 0.8, 2.4))
  expect_equal(mrpheus:::.peaks_prominence(c(1, 1, 1), 2)$peaks, integer(0))
})

test_that(".lin_interp0 equals np.interp with zero outside the range", {
  expect_equal(mrpheus:::.lin_interp0(c(0, 1, 3), c(2, 4, 0), c(-1, 0, 0.5, 1, 2.5, 3, 4)),
               c(0, 2, 3, 4, 1, 0, 0))
})

test_that(".next_fast_len equals scipy.fftpack.next_fast_len", {
  expect_equal(vapply(c(7950000, 11522304, 123457, 100001, 3000), mrpheus:::.next_fast_len, integer(1)),
               c(7962624L, 11664000L, 124416L, 101250L, 3000L))
})

test_that(".hilbert_analytic is the analytic signal of scipy.signal.hilbert", {
  n <- 1000L
  t <- (seq_len(n) - 1L) / n
  x <- sin(2 * pi * 5 * t)                       # whole periods: no leakage
  a <- mrpheus:::.hilbert_analytic(x, n)
  expect_equal(Re(a), x, tolerance = 1e-10)
  expect_equal(Im(a), -cos(2 * pi * 5 * t), tolerance = 1e-10)
  # zero padding to a longer length keeps the length of the input
  expect_length(mrpheus:::.hilbert_analytic(x, 1024L), n)
})

test_that(".mne_auto_trans follows MNE's automatic transition bands", {
  expect_equal(mrpheus:::.mne_auto_trans(100, 1, 30), c(l = 1, h = 7.5))
  expect_equal(mrpheus:::.mne_auto_trans(256, 1, 30), c(l = 1, h = 7.5))
  expect_equal(mrpheus:::.mne_auto_trans(100, 12, 15), c(l = 3, h = 3.75))
})

test_that(".sp_crit counts a non-finite value as not met, as YASA's astype(int) does", {
  expect_equal(mrpheus:::.sp_crit(c(0.1, NA, 0.9, NaN, 0.5), 0.5), c(0, 0, 1, 0, 1))
})

test_that("a flat stretch does not remove the spindles after it", {
  # A moving correlation over a window without variance is NaN (YASA returns NaN there and counts it
  # as 0). As NA it turned the cumulative sum of the smoothing into NA, and no spindle after the
  # first such sample was returned.
  psg  <- .make_sp_psg()
  base <- suppressMessages(compute_spindles(psg))
  psg$epochs[[2]]$C3[2201:2800] <- 0            # 52-58 s, between the bursts at 40 s and 70 s
  gap  <- suppressMessages(compute_spindles(psg))
  expect_equal(nrow(gap), nrow(base))
  expect_equal(gap$epoch, base$epoch)
  expect_equal(gap$start_s, base$start_s, tolerance = 0.05)
})
