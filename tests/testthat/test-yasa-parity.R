# tests/testthat/test-yasa-parity.R
#
# stage_epochs() must reproduce YASA. Fixtures from YASA 0.8.1 on a
# deterministic synthetic recording (data-raw/generate_yasa_parity_fixture.py):
# the 149 features for EEG + EOG + EMG, and posterior probabilities for all
# eight pre-trained classifiers. The recording is generated here from the same
# formulas, so no signal files are stored.
#
# Two tests, because they catch different things:
#
# 1. FEATURES (the primary guard). Every feature must equal YASA's to float32
#    precision. A feature error is large and specific: the old Simpson rule
#    moved the 0.4-1 and 1-4 Hz band powers by up to half a standard deviation,
#    the old skew / kurtosis scaling by 7e-4 relative, in columns the
#    probabilities mostly hide.
#
# 2. PROBABILITIES (an end-to-end check of classifier choice, column order,
#    metadata and float32 handling). It cannot be tight, for this reason: at
#    the recording's first and last sample the zero-phase FIR filter's output is
#    zero by construction and its sign is floating-point noise (+6e-15 in MNE,
#    -4e-14 in R here; it varies with platform). That can change a zero-crossing
#    count by one, and because `robust_scale` uses whole-column quantiles it
#    shifts every normalised zero-crossing value slightly. On this 30-epoch
#    fixture, forcing +/-1 on the first or last count in any channel moves a
#    probability by up to 2e-2 (about 1e-3 on a 780-epoch night). That is the
#    same size as the effect of the old skew / kurtosis scaling (3e-2) and
#    Simpson errors (7e-2), so probabilities cannot separate noise from those
#    bugs; test 1 does. Zero-crossing columns are therefore the only ones
#    test 1 treats specially: the normalised ones, and the raw count in the
#    first and last epoch, are skipped.

.lcg <- function(n, seed) {          # Park-Miller; 16807 * state < 2^53, exact in doubles
  out <- numeric(n); state <- seed
  for (i in seq_len(n)) {
    state <- (16807 * state) %% 2147483647
    out[i] <- state / 2147483647 - 0.5
  }
  out
}

.parity_signals <- function(fs = 100L, n_epochs = 30L) {
  n  <- fs * 30L * n_epochs
  t  <- (seq_len(n) - 1L) / fs
  tp <- 2 * pi
  env <- 1 + 0.6 * sin(tp * t / 1200)
  list(
    EEG = 30 * env * (sin(tp * 1.1 * t + 0.3) + 0.5 * sin(tp * 10.3 * t) +
                      0.3 * sin(tp * 2.7 * t + 1)) + 12 * .lcg(n, 12345),
    EOG = 20 * (sin(tp * 0.4 * t) * (1 + 0.8 * sin(tp * t / 900)) +
                0.3 * sin(tp * 5.1 * t)) + 8 * .lcg(n, 67890),
    EMG = 4 * (1 + 0.9 * sin(tp * t / 700 + 2)) * 2 * .lcg(n, 424242) +
          2 * sin(tp * 13.7 * t)
  )
}

.parity_psg <- function() {
  sig <- .parity_signals()
  edf <- structure(
    list(
      header     = list(),
      signals    = lapply(sig, function(s) list(signal = s)),
      channels   = tibble::tibble(label = names(sig), sample_rate = 100),
      duration_s = length(sig$EEG) / 100,
      path       = NA_character_
    ),
    class = "mrpheus_edf"
  )
  suppressMessages(prepare_psg(edf))
}

test_that("the synthetic recording matches the one the fixture was made from", {
  sig <- .parity_signals()
  # printed by data-raw/generate_yasa_parity_fixture.py
  expect_equal(sig$EEG[1:2], c(11.598243405202565, 32.16198145694556), tolerance = 1e-10)
  expect_equal(sig$EOG[1:2], c(0.25065767217923796, -0.8039545861659256), tolerance = 1e-10)
  expect_equal(sig$EMG[1:2], c(-2.614455640828706, 6.670863492590634), tolerance = 1e-10)
})

test_that("mrpheus features equal YASA's features to float32 precision", {
  yf  <- utils::read.csv(testthat::test_path("fixtures", "yasa_parity_features.csv"),
                         check.names = FALSE, row.names = 1)
  psg <- .parity_psg()
  mf  <- as.matrix(mrpheus:::.extract_staging_features(psg, "EEG", "EOG", "EMG")[, -1])
  y   <- as.matrix(yf)

  expect_setequal(colnames(mf), colnames(y))
  expect_equal(nrow(mf), nrow(y))
  mf <- mf[, colnames(y)]
  mf[] <- mrpheus:::.as_float32(mf)                 # YASA stores float32
  expect_equal(sum(is.na(mf)), sum(is.na(y)))

  bad <- abs(mf - y) > 1e-5 * pmax(1, abs(y))       # float32 is ~1e-7 relative; the bugs were >= 7e-4
  bad[is.na(bad)] <- FALSE
  bad[, grepl("_nzc_(c7min|p2min)_norm$", colnames(y))] <- FALSE   # see the note at the top
  bad[c(1L, nrow(y)), grepl("_nzc$", colnames(y))]      <- FALSE
  expect_false(any(bad),
               info = paste("differing:", paste(colnames(y)[colSums(bad) > 0], collapse = ", ")))
})

test_that("stage_epochs gives YASA's stage and probabilities for all eight classifiers", {
  fix <- utils::read.csv(testthat::test_path("fixtures", "yasa_parity_probabilities.csv"))
  psg <- .parity_psg()
  pcols  <- paste0("prob_", c("N1", "N2", "N3", "REM", "W"))
  labels <- c("N1", "N2", "N3", "REM", "W")

  for (key in unique(fix$classifier)) {
    use_eog <- grepl("eog", key); use_emg <- grepl("emg", key); demo <- grepl("demo", key)
    st <- suppressMessages(stage_epochs(
      psg, eeg_channel = "EEG",
      eog_channel = if (use_eog) "EOG" else NA_character_,
      emg_channel = if (use_emg) "EMG" else NA_character_,
      age  = if (demo) 50 else NULL,
      male = if (demo) 1  else NULL
    ))
    ref <- fix[fix$classifier == key, ]
    expect_equal(nrow(st), nrow(ref), info = key)
    expect_equal(attr(st, "classifier"), paste0("clf_", gsub("_", "+", key)), info = key)

    m <- as.matrix(st[, pcols]); y <- as.matrix(ref[, pcols])
    # Coarse on purpose (see the note at the top): wrong classifier, column
    # order or metadata handling moves probabilities by far more than this.
    expect_lt(max(abs(m - y)), 5e-2, label = paste(key, "max |dP|"))

    top2  <- t(apply(y, 1, function(r) sort(r, decreasing = TRUE)[1:2]))
    clear <- (top2[, 1] - top2[, 2]) > 0.1
    expect_identical(st$stage[clear], labels[apply(y, 1, which.max)][clear], info = key)
  }
})

test_that("an explicit model_path gives the same result as the automatic choice", {
  psg  <- .parity_psg()
  auto <- suppressMessages(stage_epochs(psg, "EEG", "EOG", "EMG"))
  path <- system.file("models", "yasa_eeg_eog_emg.txt.gz", package = "mrpheus")
  expl <- suppressMessages(stage_epochs(psg, "EEG", "EOG", "EMG", model_path = path))
  expect_equal(as.data.frame(auto), as.data.frame(expl), ignore_attr = TRUE)
})

test_that("age and male must be supplied together and be valid", {
  psg <- .parity_psg()
  expect_error(suppressMessages(stage_epochs(psg, "EEG", "EOG", "EMG", age = 40)), "both")
  expect_error(suppressMessages(stage_epochs(psg, "EEG", "EOG", "EMG", male = 1)), "both")
  expect_error(suppressMessages(stage_epochs(psg, "EEG", "EOG", "EMG", age = 150, male = 1)), "between 0 and 120")
  expect_error(suppressMessages(stage_epochs(psg, "EEG", "EOG", "EMG", age = 40, male = 2)), "0 or 1")
})

test_that("the classifier key follows YASA's naming", {
  expect_equal(mrpheus:::.yasa_model_key(FALSE, FALSE, FALSE), "eeg")
  expect_equal(mrpheus:::.yasa_model_key(TRUE,  FALSE, FALSE), "eeg_eog")
  expect_equal(mrpheus:::.yasa_model_key(FALSE, TRUE,  FALSE), "eeg_emg")
  expect_equal(mrpheus:::.yasa_model_key(TRUE,  TRUE,  TRUE),  "eeg_eog_emg_demo")
  for (key in c("eeg", "eeg_eog", "eeg_emg", "eeg_eog_emg", "eeg_demo", "eeg_eog_demo",
                "eeg_emg_demo", "eeg_eog_emg_demo")) {
    expect_true(nzchar(system.file("models", paste0("yasa_", key, ".txt.gz"), package = "mrpheus")),
                info = key)
  }
})

test_that(".as_float32 rounds like numpy's astype(float32)", {
  f32 <- mrpheus:::.as_float32
  expect_identical(f32(0.1), 0.10000000149011612)       # float32(0.1) as a double
  expect_identical(f32(16777217), 16777216)             # 2^24 + 1 is not representable
  expect_identical(f32(c(1, 2.5, -3)), c(1, 2.5, -3))   # exactly representable values are unchanged
  expect_true(is.nan(f32(NA_real_)))
  expect_equal(dim(matrix(f32(matrix(0.1, 2, 2)), 2, 2)), c(2L, 2L))
})
