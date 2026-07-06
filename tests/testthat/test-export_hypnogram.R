# tests/testthat/test-export_hypnogram.R
#
# Tests for export_hypnogram() and .sleep_architecture().
# Expected values validated against yasa.sleep_statistics() directly.

# ── Reference hypnogram (from YASA docstring) ─────────────────────────────────
# hypno = [0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 2, 3, 3, 4, 4, 4, 4, 0, 0]
# sf_hyp = 1/30  =>  epoch_s = 30
# Expected output (YASA):
#   TIB=10, SPT=8, WASO=0, TST=8, N1=1.5, N2=2, N3=2.5, REM=2, NREM=6
#   SOL=1, Lat_N1=1, Lat_N2=2.5, Lat_N3=4, Lat_REM=7
#   %N1=18.75, %N2=25, %N3=31.25, %REM=25, %NREM=75
#   SE=80, SME=100

YASA_HYPNO <- as.integer(
  c(0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 2, 3, 3, 4, 4, 4, 4, 0, 0)
)

make_staging <- function(hypno = YASA_HYPNO) {
  tibble::tibble(epoch = seq_along(hypno), stage = hypno)
}

# ── .sleep_architecture parity with YASA ─────────────────────────────────────

test_that(".sleep_architecture matches YASA reference values", {
  arch <- mrpheus:::.sleep_architecture(YASA_HYPNO, epoch_s = 30)

  expect_equal(arch$TIB,  10.0,   tolerance = 1e-10)
  expect_equal(arch$SPT,   8.0,   tolerance = 1e-10)
  expect_equal(arch$WASO,  0.0,   tolerance = 1e-10)
  expect_equal(arch$TST,   8.0,   tolerance = 1e-10)
  expect_equal(arch$N1,    1.5,   tolerance = 1e-10)
  expect_equal(arch$N2,    2.0,   tolerance = 1e-10)
  expect_equal(arch$N3,    2.5,   tolerance = 1e-10)
  expect_equal(arch$REM,   2.0,   tolerance = 1e-10)
  expect_equal(arch$NREM,  6.0,   tolerance = 1e-10)
  expect_equal(arch$SOL,   1.0,   tolerance = 1e-10)
  expect_equal(arch$Lat_N1,  1.0, tolerance = 1e-10)
  expect_equal(arch$Lat_N2,  2.5, tolerance = 1e-10)
  expect_equal(arch$Lat_N3,  4.0, tolerance = 1e-10)
  expect_equal(arch$Lat_REM, 7.0, tolerance = 1e-10)
  expect_equal(arch$pct_N1,  18.75, tolerance = 1e-10)
  expect_equal(arch$pct_N2,  25.0,  tolerance = 1e-10)
  expect_equal(arch$pct_N3,  31.25, tolerance = 1e-10)
  expect_equal(arch$pct_REM, 25.0,  tolerance = 1e-10)
  expect_equal(arch$pct_NREM, 75.0, tolerance = 1e-10)
  expect_equal(arch$SE,  80.0,  tolerance = 1e-10)
  expect_equal(arch$SME, 100.0, tolerance = 1e-10)
})

test_that(".sleep_architecture Lat_REM_AASM = Lat_REM - SOL", {
  arch <- mrpheus:::.sleep_architecture(YASA_HYPNO, epoch_s = 30)
  expect_equal(arch$Lat_REM_AASM, arch$Lat_REM - arch$SOL, tolerance = 1e-10)
})

test_that(".sleep_architecture handles all-wake hypnogram", {
  hypno <- as.integer(rep(0L, 20L))
  arch  <- mrpheus:::.sleep_architecture(hypno, epoch_s = 30)
  expect_equal(arch$TST, 0)
  expect_equal(arch$SE,  0)
  expect_true(is.na(arch$SOL))
  expect_true(is.na(arch$WASO))
  expect_true(is.na(arch$Lat_REM))
})

test_that(".sleep_architecture handles no REM", {
  hypno <- as.integer(c(0, 0, 1, 2, 2, 3, 3, 0))
  arch  <- mrpheus:::.sleep_architecture(hypno, epoch_s = 30)
  expect_true(is.na(arch$Lat_REM))
  expect_true(is.na(arch$Lat_REM_AASM))
  # pct_REM = 100 * 0 / TST = 0 (matches YASA behaviour)
  expect_equal(arch$pct_REM, 0, tolerance = 1e-10)
})

test_that(".sleep_architecture excludes artefact epochs from TST", {
  # Artefact (-1) between sleep stages should not count toward TST
  hypno <- as.integer(c(0, 1, -1, 2, 3, 4))
  arch  <- mrpheus:::.sleep_architecture(hypno, epoch_s = 30)
  # TST = sleep epochs within SPT = epochs >0 only: 1,2,3,4 = 4 epochs
  expect_equal(arch$TST, 4 * 30 / 60, tolerance = 1e-10)
})

test_that(".sleep_architecture WASO counts wake within SPT only", {
  # Wake at start and end is not WASO
  hypno <- as.integer(c(0, 0, 2, 0, 2, 0, 0))
  arch  <- mrpheus:::.sleep_architecture(hypno, epoch_s = 30)
  # SPT = epochs 3–5, WASO = 1 wake epoch within SPT
  expect_equal(arch$WASO, 1 * 30 / 60, tolerance = 1e-10)
})

# ── export_hypnogram ──────────────────────────────────────────────────────────

test_that("export_hypnogram returns mrpheus_hypnogram class", {
  out <- export_hypnogram(make_staging())
  expect_s3_class(out, "mrpheus_hypnogram")
})

test_that("export_hypnogram attaches sleep_architecture attribute", {
  out  <- export_hypnogram(make_staging())
  arch <- attr(out, "sleep_architecture")
  expect_type(arch, "list")
  expect_true(all(c("TIB", "TST", "SE", "SOL", "Lat_REM_AASM") %in% names(arch)))
})

test_that("export_hypnogram preserves original columns", {
  staging <- make_staging()
  out     <- export_hypnogram(staging)
  expect_true(all(names(staging) %in% names(out)))
})

test_that("export_hypnogram stores metadata attributes correctly", {
  out <- export_hypnogram(
    make_staging(),
    epoch_s        = 30,
    participant_id = "sub-001",
    start_time     = as.POSIXct("2024-01-01 22:00:00")
  )
  expect_equal(attr(out, "epoch_s"),        30)
  expect_equal(attr(out, "participant_id"), "sub-001")
  expect_equal(attr(out, "source"),         "mrpheus")
  expect_equal(attr(out, "resolution"),     "AASM")
  expect_s3_class(attr(out, "start_time"), "POSIXct")
})

test_that("export_hypnogram errors on missing stage column", {
  expect_error(
    export_hypnogram(tibble::tibble(epoch = 1:5)),
    regexp = "stage"
  )
})
