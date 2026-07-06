# tests/testthat/test-correct_eog.R
#
# Tests for correct_eog_regression() and correct_eog_ica().
#
# Approach: synthetic PSG with a known EOG artifact added to the EEG.
# After correction the EEG-EOG correlation must be lower than before.

# ── Helper ───────────────────────────────────────────────────────────────────

make_eog_psg <- function(sr = 100L, n_epochs = 6L, alpha = 0.8, seed = 42L) {
  set.seed(seed)
  epoch_s  <- 30L
  n_ep     <- as.integer(sr * epoch_s)
  n_tot    <- n_ep * n_epochs
  t        <- seq_len(n_tot) / sr

  # EOG: blink-like signal (low-frequency sinusoid + spikes)
  eog <- sin(2 * pi * 0.3 * t) + 0.5 * sin(2 * pi * 0.8 * t)

  # Two EEG channels: broadband neural + EOG contamination
  neural1 <- rnorm(n_tot) * 0.5
  neural2 <- rnorm(n_tot) * 0.5
  eeg1    <- neural1 + alpha * eog
  eeg2    <- neural2 + alpha * eog

  edf <- list(
    channels = data.frame(
      label       = c("EEG C3", "EEG C4", "EOG LOC"),
      sample_rate = c(sr, sr, sr),
      stringsAsFactors = FALSE
    ),
    signals = list(
      `EEG C3`  = list(signal = eeg1),
      `EEG C4`  = list(signal = eeg2),
      `EOG LOC` = list(signal = eog)
    )
  )
  cmap <- tibble::tibble(
    label       = c("EEG C3", "EEG C4", "EOG LOC"),
    type        = c("EEG", "EEG", "EOG"),
    sample_rate = c(sr, sr, sr),
    bad         = c(FALSE, FALSE, FALSE)
  )
  ep_n   <- as.integer(epoch_s * sr)
  epochs <- lapply(seq_len(n_epochs), function(i) {
    start <- (i - 1L) * ep_n + 1L
    list(
      `EEG C3`  = eeg1[start:(start + ep_n - 1L)],
      `EEG C4`  = eeg2[start:(start + ep_n - 1L)],
      `EOG LOC` = eog[start:(start + ep_n - 1L)]
    )
  })
  structure(
    list(edf = edf, epochs = epochs, n_epochs = n_epochs,
         epoch_s = epoch_s, channel_map = cmap),
    class = "mrpheus_psg"
  )
}

# ── correct_eog_regression ────────────────────────────────────────────────────

test_that("correct_eog_regression returns an mrpheus_psg", {
  psg <- make_eog_psg()
  out <- correct_eog_regression(psg, verbose = FALSE)
  expect_s3_class(out, "mrpheus_psg")
})

test_that("correct_eog_regression preserves n_epochs and epoch_s", {
  psg <- make_eog_psg()
  out <- correct_eog_regression(psg, verbose = FALSE)
  expect_equal(out$n_epochs, psg$n_epochs)
  expect_equal(out$epoch_s,  psg$epoch_s)
})

test_that("correct_eog_regression preserves epoch signal lengths", {
  psg <- make_eog_psg()
  out <- correct_eog_regression(psg, verbose = FALSE)
  for (i in seq_len(out$n_epochs)) {
    expect_equal(length(out$epochs[[i]][["EEG C3"]]),
                 length(psg$epochs[[i]][["EEG C3"]]))
  }
})

test_that("correct_eog_regression reduces EEG-EOG correlation", {
  psg <- make_eog_psg(alpha = 0.9)
  eeg_before <- psg$edf$signals[["EEG C3"]]$signal
  eog        <- psg$edf$signals[["EOG LOC"]]$signal

  out        <- correct_eog_regression(psg, verbose = FALSE)
  eeg_after  <- out$edf$signals[["EEG C3"]]$signal

  r_before <- abs(cor(eeg_before, eog))
  r_after  <- abs(cor(eeg_after,  eog))
  expect_lt(r_after, r_before)
})

test_that("correct_eog_regression leaves EOG channel unchanged", {
  psg <- make_eog_psg()
  out <- correct_eog_regression(psg, verbose = FALSE)
  expect_equal(out$edf$signals[["EOG LOC"]]$signal,
               psg$edf$signals[["EOG LOC"]]$signal)
})

test_that("correct_eog_regression errors when no EOG channels found", {
  psg  <- make_eog_psg()
  # Mark EOG as EEG so auto-detection finds no EOG
  psg$channel_map$type[psg$channel_map$type == "EOG"] <- "EEG"
  expect_error(correct_eog_regression(psg, verbose = FALSE), regexp = "EOG")
})

test_that("correct_eog_regression errors on non-mrpheus_psg input", {
  expect_error(correct_eog_regression(list()), regexp = "mrpheus_psg")
})

# ── correct_eog_ica ───────────────────────────────────────────────────────────

test_that("correct_eog_ica returns an mrpheus_psg", {
  psg <- make_eog_psg()
  set.seed(1)
  out <- correct_eog_ica(psg, verbose = FALSE)
  expect_s3_class(out, "mrpheus_psg")
})

test_that("correct_eog_ica preserves n_epochs and epoch_s", {
  psg <- make_eog_psg()
  set.seed(1)
  out <- correct_eog_ica(psg, verbose = FALSE)
  expect_equal(out$n_epochs, psg$n_epochs)
  expect_equal(out$epoch_s,  psg$epoch_s)
})

test_that("correct_eog_ica preserves epoch signal lengths", {
  psg <- make_eog_psg()
  set.seed(1)
  out <- correct_eog_ica(psg, verbose = FALSE)
  for (i in seq_len(out$n_epochs)) {
    expect_equal(length(out$epochs[[i]][["EEG C3"]]),
                 length(psg$epochs[[i]][["EEG C3"]]))
  }
})

test_that("correct_eog_ica reduces EEG-EOG correlation", {
  psg <- make_eog_psg(alpha = 0.95, n_epochs = 10L, seed = 7L)
  eeg_before <- psg$edf$signals[["EEG C3"]]$signal
  eog        <- psg$edf$signals[["EOG LOC"]]$signal

  set.seed(42)
  out       <- correct_eog_ica(psg, threshold = 0.3, verbose = FALSE)
  eeg_after <- out$edf$signals[["EEG C3"]]$signal

  r_before <- abs(cor(eeg_before, eog))
  r_after  <- abs(cor(eeg_after,  eog))
  expect_lt(r_after, r_before)
})

test_that("correct_eog_ica leaves EOG channel unchanged", {
  psg <- make_eog_psg()
  set.seed(1)
  out <- correct_eog_ica(psg, verbose = FALSE)
  expect_equal(out$edf$signals[["EOG LOC"]]$signal,
               psg$edf$signals[["EOG LOC"]]$signal)
})

test_that("correct_eog_ica errors on single EEG channel", {
  psg <- make_eog_psg()
  expect_error(
    correct_eog_ica(psg,
                    eeg_channels = "EEG C3",  # only 1 channel
                    verbose      = FALSE),
    regexp = "2 EEG"
  )
})

test_that("correct_eog_ica returns unchanged psg when threshold not exceeded", {
  psg <- make_eog_psg(alpha = 0.01)  # near-zero EOG contamination
  set.seed(1)
  # Very high threshold: no component should be flagged
  out <- correct_eog_ica(psg, threshold = 0.99, verbose = FALSE)
  expect_equal(out$edf$signals[["EEG C3"]]$signal,
               psg$edf$signals[["EEG C3"]]$signal)
})

test_that("correct_eog_ica errors when no EOG channels found", {
  psg <- make_eog_psg()
  psg$channel_map$type[psg$channel_map$type == "EOG"] <- "EEG"
  set.seed(1)
  expect_error(correct_eog_ica(psg, verbose = FALSE), regexp = "EOG")
})

test_that("correct_eog_ica errors on non-mrpheus_psg input", {
  expect_error(correct_eog_ica(list(), rep(0L, 5L)), regexp = "mrpheus_psg")
})
