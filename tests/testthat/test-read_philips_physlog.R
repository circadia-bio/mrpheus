fixture_path <- function(name) {
  testthat::test_path("fixtures", name)
}

fixture_exists <- function(name) {
  file.exists(fixture_path(name))
}

read_fixture <- function(name) {
  utils::read.csv(fixture_path(name), header = FALSE)
}

example_physlog <- function(...) {
  path <- system.file("extdata", "example_physlog.log", package = "mrpheus")
  read_philips_physlog(path, ...)
}

# ---------------------------------------------------------------------------

test_that("read_philips_physlog returns correct S3 class and structure", {
  rec <- example_physlog()

  expect_s3_class(rec, "mrpheus_physlog")
  expect_named(rec, c("C", "M", "I", "HDR", "path"))
})

test_that("read_philips_physlog signal matrix has correct dimensions and type", {
  rec <- example_physlog()

  expect_true(is.integer(rec$C))
  expect_equal(ncol(rec$C), 10L)
  expect_equal(nrow(rec$C), 496L * 20L)
  expect_equal(
    colnames(rec$C),
    c("v1raw", "v2raw", "v1", "v2", "ppu", "resp", "gx", "gy", "gz", "mark")
  )
})

test_that("read_philips_physlog HDR carries resolved sfreq and system", {
  rec_wBTU   <- example_physlog(system = "wBTU")
  rec_wired  <- example_physlog(system = "wired")
  rec_custom <- example_physlog(system = "custom", sfreq = 400)

  expect_equal(rec_wBTU$HDR$sfreq,   496)
  expect_equal(rec_wBTU$HDR$system,  "wBTU")
  expect_equal(rec_wired$HDR$sfreq,  500)
  expect_equal(rec_custom$HDR$sfreq, 400)
})

test_that("read_philips_physlog extracts scanner markers at correct positions", {
  rec <- example_physlog()
  n   <- nrow(rec$C)

  expect_equal(rec$I$ScannerStart, 496L)
  expect_equal(rec$I$ScannerStop,  n - 496L + 42L)
})

test_that("read_philips_physlog extracts VcgOnset markers", {
  rec <- example_physlog()

  expect_equal(length(rec$I$VcgOnset), 50L)
  expect_true(all(rec$I$VcgOnset >= 1L))
  expect_true(all(rec$I$VcgOnset <= nrow(rec$C)))
})

test_that("read_philips_physlog channel selection works", {
  path <- system.file("extdata", "example_physlog.log", package = "mrpheus")

  rec_ecg  <- read_philips_physlog(path, channels = c("v1raw", "v2raw"))
  expect_equal(ncol(rec_ecg$C),     2L)
  expect_equal(colnames(rec_ecg$C), c("v1raw", "v2raw"))

  rec_none <- read_philips_physlog(path, channels = "none")
  expect_equal(ncol(rec_none$C), 0L)
})

test_that("read_philips_physlog errors on missing file", {
  expect_error(read_philips_physlog("nonexistent.log"), regexp = "not found")
})

test_that("read_philips_physlog errors on missing sfreq for custom system", {
  path <- system.file("extdata", "example_physlog.log", package = "mrpheus")
  expect_error(read_philips_physlog(path, system = "custom"), regexp = "sfreq")
})

# ---------------------------------------------------------------------------
# MATLAB bit-perfect comparison
# Skipped until fixtures are generated via data-raw/generate_matlab_fixtures.m
# ---------------------------------------------------------------------------

test_that("signal matrix matches MATLAB ReadPhilipsScanPhysLog — first 100 rows", {
  skip_if(
    !fixture_exists("fixture_C_first100.csv"),
    "MATLAB fixtures not present — run: matlab -batch \"run data-raw/generate_matlab_fixtures.m\""
  )

  rec <- example_physlog()
  ref <- as.matrix(read_fixture("fixture_C_first100.csv"))
  storage.mode(ref) <- "integer"

  expect_identical(unname(rec$C[1:100, ]), unname(ref))
})

test_that("signal matrix matches MATLAB ReadPhilipsScanPhysLog — last 100 rows", {
  skip_if(
    !fixture_exists("fixture_C_last100.csv"),
    "MATLAB fixtures not present — run: matlab -batch \"run data-raw/generate_matlab_fixtures.m\""
  )

  rec <- example_physlog()
  n   <- nrow(rec$C)
  ref <- as.matrix(read_fixture("fixture_C_last100.csv"))
  storage.mode(ref) <- "integer"

  expect_identical(unname(rec$C[(n - 99):n, ]), unname(ref))
})

test_that("markers match MATLAB ReadPhilipsScanPhysLog", {
  skip_if(
    !fixture_exists("fixture_scanner_start.csv"),
    "MATLAB fixtures not present — run: matlab -batch \"run data-raw/generate_matlab_fixtures.m\""
  )

  rec       <- example_physlog()
  ref_start <- as.integer(read_fixture("fixture_scanner_start.csv")[[1]])
  ref_stop  <- as.integer(read_fixture("fixture_scanner_stop.csv")[[1]])
  ref_vcg   <- as.integer(read_fixture("fixture_vcg_onset.csv")[[1]])

  expect_identical(rec$I$ScannerStart, ref_start)
  expect_identical(rec$I$ScannerStop,  ref_stop)
  expect_identical(rec$I$VcgOnset,     ref_vcg)
})
