fixture_path <- function(name) {
  testthat::test_path("fixtures", name)
}

fixture_exists <- function(name) {
  file.exists(fixture_path(name))
}

read_fixture <- function(name) {
  utils::read.csv(fixture_path(name), header = FALSE)
}

# ---------------------------------------------------------------------------

test_that("detect_qrs returns correct S3 class and structure", {
  path <- system.file("extdata", "example_physlog.log", package = "mrpheus")
  rec  <- read_philips_physlog(path)
  qrs  <- detect_qrs(as.double(rec$C[, "v1raw"]), fs = rec$HDR$sfreq)

  expect_s3_class(qrs, "mrpheus_qrs")
  expect_named(qrs, c("qrs_i", "qrs_amp", "delay", "fs"))
  expect_equal(qrs$fs, rec$HDR$sfreq)
  expect_true(qrs$delay > 0)
})

test_that("detect_qrs output vectors have consistent lengths and types", {
  path <- system.file("extdata", "example_physlog.log", package = "mrpheus")
  rec  <- read_philips_physlog(path)
  qrs  <- detect_qrs(as.double(rec$C[, "v1raw"]), fs = rec$HDR$sfreq)

  expect_equal(length(qrs$qrs_i), length(qrs$qrs_amp))
  expect_true(is.integer(qrs$qrs_i))
  expect_true(is.double(qrs$qrs_amp))
})

test_that("detect_qrs detects the right number of peaks in synthetic data", {
  path <- system.file("extdata", "example_physlog.log", package = "mrpheus")
  rec  <- read_philips_physlog(path)
  qrs  <- detect_qrs(as.double(rec$C[, "v1raw"]), fs = rec$HDR$sfreq)

  # Synthetic data: 50 R-peaks at 150 bpm
  expect_gte(length(qrs$qrs_i), 45L)
  expect_lte(length(qrs$qrs_i), 55L)
})

test_that("detect_qrs peak indices are within signal bounds", {
  path <- system.file("extdata", "example_physlog.log", package = "mrpheus")
  rec  <- read_philips_physlog(path)
  ecg  <- as.double(rec$C[, "v1raw"])
  qrs  <- detect_qrs(ecg, fs = rec$HDR$sfreq)

  expect_true(all(qrs$qrs_i >= 1L))
  expect_true(all(qrs$qrs_i <= length(ecg)))
})

test_that("detect_qrs produces plausible heart rate in synthetic data", {
  path    <- system.file("extdata", "example_physlog.log", package = "mrpheus")
  rec     <- read_philips_physlog(path)
  qrs     <- detect_qrs(as.double(rec$C[, "v1raw"]), fs = rec$HDR$sfreq)
  mean_hr <- 60 / (mean(diff(qrs$qrs_i)) / rec$HDR$sfreq)

  # Synthetic HR is 150 bpm; allow ±30 bpm
  expect_gte(mean_hr, 120)
  expect_lte(mean_hr, 180)
})

test_that("detect_qrs peaks align with VcgOnset markers in synthetic data", {
  path        <- system.file("extdata", "example_physlog.log", package = "mrpheus")
  rec         <- read_philips_physlog(path)
  qrs         <- detect_qrs(as.double(rec$C[, "v1raw"]), fs = rec$HDR$sfreq)
  known_peaks <- rec$I$VcgOnset

  # Each detected peak should be within 15 samples of a known marker
  matched <- vapply(
    qrs$qrs_i,
    function(p) any(abs(known_peaks - p) <= 15L),
    logical(1)
  )
  expect_gte(mean(matched), 0.9)
})

test_that("detect_qrs validates inputs", {
  expect_error(detect_qrs("not numeric", fs = 496))
  expect_error(detect_qrs(c(1, 2, 3), fs = -1))
  expect_error(detect_qrs(c(1, 2, 3), fs = c(496, 496)))
})

# ---------------------------------------------------------------------------
# MATLAB comparison — same peaks within 5-sample tolerance
# Skipped until fixtures are generated via data-raw/generate_matlab_fixtures.m
# ---------------------------------------------------------------------------

test_that("detect_qrs R-peak count matches MATLAB pan_tompkin", {
  skip_if(
    !fixture_exists("fixture_qrs_i_raw.csv"),
    "MATLAB fixtures not present — run: matlab -batch \"run data-raw/generate_matlab_fixtures.m\""
  )

  path    <- system.file("extdata", "example_physlog.log", package = "mrpheus")
  rec     <- read_philips_physlog(path)
  qrs     <- detect_qrs(as.double(rec$C[, "v1raw"]), fs = rec$HDR$sfreq)
  ref_i   <- as.integer(read_fixture("fixture_qrs_i_raw.csv")[[1]])

  expect_equal(length(qrs$qrs_i), length(ref_i))
})

test_that("detect_qrs R-peak positions match MATLAB pan_tompkin within 5 samples", {
  skip_if(
    !fixture_exists("fixture_qrs_i_raw.csv"),
    "MATLAB fixtures not present — run: matlab -batch \"run data-raw/generate_matlab_fixtures.m\""
  )

  path  <- system.file("extdata", "example_physlog.log", package = "mrpheus")
  rec   <- read_philips_physlog(path)
  qrs   <- detect_qrs(as.double(rec$C[, "v1raw"]), fs = rec$HDR$sfreq)
  ref_i <- as.integer(read_fixture("fixture_qrs_i_raw.csv")[[1]])

  diffs <- abs(sort(qrs$qrs_i) - sort(ref_i))
  expect_true(
    all(diffs <= 5L),
    info = sprintf(
      "Max deviation: %d samples; mean: %.1f samples",
      max(diffs), mean(diffs)
    )
  )
})
