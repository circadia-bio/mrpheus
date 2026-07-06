# tests/testthat/test-compute_temporal_bandpower.R
#
# Unit tests for compute_temporal_bandpower().

# ── Helper ───────────────────────────────────────────────────────────────────

make_tbp_psg <- function(sr = 100L, n_epochs = 20L, seed = 7L) {
  set.seed(seed)
  epoch_s <- 30L
  n_ep    <- as.integer(sr * epoch_s)
  n_tot   <- n_ep * n_epochs
  sig     <- rnorm(n_tot)

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
  psg <- structure(
    list(edf = edf, epochs = epochs, n_epochs = n_epochs,
         epoch_s = epoch_s, channel_map = cmap),
    class = "mrpheus_psg"
  )
  hypno <- as.integer(rep(c(0L, 2L, 3L, 4L), length.out = n_epochs))
  list(psg = psg, hypno = hypno)
}

# ── Structure ─────────────────────────────────────────────────────────────────

test_that("compute_temporal_bandpower returns a tibble", {
  d   <- make_tbp_psg()
  out <- compute_temporal_bandpower(d$psg, d$hypno)
  expect_s3_class(out, "tbl_df")
})

test_that("compute_temporal_bandpower has the expected columns", {
  d   <- make_tbp_psg()
  out <- compute_temporal_bandpower(d$psg, d$hypno)
  expect_true(all(c("time_hours", "epoch_start", "epoch_end",
                    "dominant_stage", "band", "power",
                    "relative_power") %in% names(out)))
})

test_that("compute_temporal_bandpower returns one row per window per band", {
  d           <- make_tbp_psg(n_epochs = 20L)
  n_bands     <- 6L
  window_ep   <- 10L
  step_ep     <- 5L
  n_windows   <- length(seq(1L, 20L - window_ep + 1L, by = step_ep))
  out <- compute_temporal_bandpower(d$psg, d$hypno,
                                     window_epochs = window_ep,
                                     step_epochs   = step_ep)
  expect_equal(nrow(out), n_windows * n_bands)
})

# ── Timing ───────────────────────────────────────────────────────────────────

test_that("time_hours starts at 0 and advances correctly", {
  d   <- make_tbp_psg(n_epochs = 20L)
  out <- compute_temporal_bandpower(d$psg, d$hypno,
                                     window_epochs = 10L, step_epochs = 5L)
  times <- unique(out$time_hours)
  expect_equal(times[1], 0)
  # Second window starts at epoch 6 (0-indexed: epoch 5): 5 * 30 / 3600
  expect_equal(times[2], 5 * 30 / 3600, tolerance = 1e-10)
})

test_that("epoch_start and epoch_end are within valid range", {
  d   <- make_tbp_psg(n_epochs = 20L)
  out <- compute_temporal_bandpower(d$psg, d$hypno)
  expect_true(all(out$epoch_start >= 1L))
  expect_true(all(out$epoch_end   <= d$psg$n_epochs))
  expect_true(all(out$epoch_start <= out$epoch_end))
})

# ── Dominant stage ────────────────────────────────────────────────────────────

test_that("dominant_stage is a valid stage code or NA", {
  d   <- make_tbp_psg()
  out <- compute_temporal_bandpower(d$psg, d$hypno)
  valid <- out$dominant_stage %in% c(0L, 1L, 2L, 3L, 4L, NA_integer_)
  expect_true(all(valid))
})

test_that("dominant_stage is NA when all epochs are artefacts", {
  d         <- make_tbp_psg(n_epochs = 20L)
  hypno_art <- rep(-1L, 20L)
  out <- compute_temporal_bandpower(d$psg, hypno_art)
  expect_true(all(is.na(out$dominant_stage)))
})

test_that("dominant_stage reflects the most common stage", {
  d         <- make_tbp_psg(n_epochs = 20L)
  # First 10 epochs: 8 x N2 (2), 2 x Wake (0) -> dominant = 2
  hypno_fixed         <- rep(0L, 20L)
  hypno_fixed[1:8]    <- 2L
  out <- compute_temporal_bandpower(d$psg, hypno_fixed,
                                     window_epochs = 10L, step_epochs = 10L)
  expect_equal(out$dominant_stage[1], 2L)
})

# ── Band power values ─────────────────────────────────────────────────────────

test_that("power values are non-negative", {
  d   <- make_tbp_psg()
  out <- compute_temporal_bandpower(d$psg, d$hypno)
  expect_true(all(out$power >= 0, na.rm = TRUE))
})

test_that("relative_power sums to 1 within each time window", {
  d   <- make_tbp_psg()
  out <- compute_temporal_bandpower(d$psg, d$hypno)
  totals <- tapply(out$relative_power, out$time_hours, sum, na.rm = TRUE)
  expect_equal(as.vector(totals), rep(1, length(totals)), tolerance = 1e-10)
})

test_that("delta band dominates a delta-enriched signal", {
  sr      <- 100L
  n_ep    <- 20L
  epoch_s <- 30L
  n_samp  <- as.integer(sr * epoch_s)
  set.seed(99)
  t <- seq_len(n_samp * n_ep) / sr
  sig <- sin(2 * pi * 2 * t) + 0.01 * rnorm(length(t))  # 2 Hz dominates

  edf <- list(
    channels = data.frame(label = "EEG C3", sample_rate = sr,
                          stringsAsFactors = FALSE),
    signals  = list(`EEG C3` = list(signal = sig))
  )
  cmap <- tibble::tibble(label = "EEG C3", type = "EEG",
                          sample_rate = sr, bad = FALSE)
  epochs <- lapply(seq_len(n_ep), function(i) {
    start <- (i - 1L) * n_samp + 1L
    list(`EEG C3` = sig[start:(start + n_samp - 1L)])
  })
  psg_d <- structure(
    list(edf = edf, epochs = epochs, n_epochs = n_ep,
         epoch_s = epoch_s, channel_map = cmap),
    class = "mrpheus_psg"
  )
  hypno <- rep(3L, n_ep)  # all N3
  out <- compute_temporal_bandpower(psg_d, hypno,
                                     window_epochs = 10L, step_epochs = 10L)
  # Delta should have highest relative power in each window
  for (t_val in unique(out$time_hours)) {
    win <- out[out$time_hours == t_val, ]
    expect_equal(win$band[which.max(win$relative_power)], "delta")
  }
})

# ── Custom parameters ─────────────────────────────────────────────────────────

test_that("custom bands are respected", {
  d   <- make_tbp_psg()
  out <- compute_temporal_bandpower(d$psg, d$hypno,
                                     bands = list(slow = c(0.5, 2),
                                                  fast = c(20, 40)))
  expect_setequal(unique(out$band), c("slow", "fast"))
})

test_that("step_epochs = window_epochs gives non-overlapping windows", {
  d         <- make_tbp_psg(n_epochs = 20L)
  out       <- compute_temporal_bandpower(d$psg, d$hypno,
                                           window_epochs = 10L,
                                           step_epochs   = 10L)
  n_windows <- length(unique(out$time_hours))
  expect_equal(n_windows, 2L)  # 20 epochs / 10 = 2 non-overlapping windows
})

# ── Error handling ────────────────────────────────────────────────────────────

test_that("mismatched hypno length raises an error", {
  d <- make_tbp_psg()
  expect_error(
    compute_temporal_bandpower(d$psg, d$hypno[-1]),
    regexp = "n_epochs"
  )
})

test_that("window_epochs larger than n_epochs raises an error", {
  d <- make_tbp_psg(n_epochs = 5L)
  expect_error(
    compute_temporal_bandpower(d$psg, d$hypno, window_epochs = 10L),
    regexp = "window_epochs"
  )
})

test_that("non-mrpheus_psg input raises an error", {
  expect_error(
    compute_temporal_bandpower(list(), rep(0L, 10L)),
    regexp = "mrpheus_psg"
  )
})
