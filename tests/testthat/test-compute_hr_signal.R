test_that("compute_hr_signal returns correct values for known input", {
  # Two R-peaks 496 samples apart at 496 Hz = 1 s interval = 60 bpm
  qrs <- structure(
    list(qrs_i = c(100L, 596L), qrs_amp = c(1, 1), delay = 37, fs = 496),
    class = "mrpheus_qrs"
  )
  hr <- compute_hr_signal(qrs)

  expect_length(hr, 596L)
  expect_true(all(hr[100:596] == 60))
  expect_true(all(hr[1:99]   == 0))
})

test_that("compute_hr_signal returns correct values for 150 bpm", {
  sfreq    <- 496
  rr_samp  <- round(60 / 150 * sfreq)   # 198 samples
  r_peaks  <- c(rr_samp, 2 * rr_samp)

  qrs <- structure(
    list(qrs_i = as.integer(r_peaks), qrs_amp = c(1, 1), delay = 37,
         fs = sfreq),
    class = "mrpheus_qrs"
  )
  hr <- compute_hr_signal(qrs)

  expected_hr <- 60 * sfreq / rr_samp
  expect_equal(unique(hr[hr > 0]), expected_hr, tolerance = 1e-9)
})

test_that("compute_hr_signal output length equals last R-peak index", {
  path <- system.file("extdata", "example_physlog.log", package = "mrpheus")
  rec  <- read_philips_physlog(path)
  qrs  <- detect_qrs(as.double(rec$C[, "v1raw"]), fs = rec$HDR$sfreq)
  hr   <- compute_hr_signal(qrs)

  expect_equal(length(hr), max(qrs$qrs_i))
})

test_that("compute_hr_signal HR values are physiologically plausible", {
  path <- system.file("extdata", "example_physlog.log", package = "mrpheus")
  rec  <- read_philips_physlog(path)
  qrs  <- detect_qrs(as.double(rec$C[, "v1raw"]), fs = rec$HDR$sfreq)
  hr   <- compute_hr_signal(qrs)

  hr_nonzero <- hr[hr > 0]
  expect_gte(min(hr_nonzero), 80)
  expect_lte(max(hr_nonzero), 250)
})

test_that("compute_hr_signal pre-beat samples are zero", {
  qrs <- structure(
    list(qrs_i = c(200L, 400L), qrs_amp = c(1, 1), delay = 37, fs = 496),
    class = "mrpheus_qrs"
  )
  hr <- compute_hr_signal(qrs)

  expect_true(all(hr[1:199] == 0))
  expect_true(all(hr[200:400] > 0))
})

test_that("compute_hr_signal warns on fewer than two peaks", {
  qrs <- structure(
    list(qrs_i = 100L, qrs_amp = 1, delay = 37, fs = 496),
    class = "mrpheus_qrs"
  )
  expect_warning(compute_hr_signal(qrs), regexp = "two R-peaks")
})

test_that("compute_hr_signal errors on bad input", {
  expect_error(compute_hr_signal(list(x = 1)))
  expect_error(compute_hr_signal("not a list"))
})
