# tests/testthat/test-yasa-bandpower-parity.R
#
# compute_band_power() with its default arguments must reproduce
# yasa.bandpower() (YASA 0.8.1, default arguments). The fixture
# (tests/testthat/fixtures/yasa_bandpower_parity.csv, made by
# data-raw/generate_yasa_bandpower_parity_fixture.py) holds absolute band power,
# relative band power and TotalAbsPow for six 30 s epochs of a deterministic
# synthetic EEG at 100 Hz and at 256 Hz. The signal is generated here from the
# same formulas, so no signal files are stored.
#
# What the check guards. Each old behaviour was put back, one at a time, and the
# tests fail (the size is the median change in YASA's own band power, measured
# on Sleep-EDF SC4001E0):
#   median averaging (was mean)            ~10 %
#   Hamming window (was Hann)              ~2 %
#   periodic window (was symmetric)        up to ~0.8 %
#   nfft = window (was 1024)               ~3 % at 100 Hz
#   half-window overlap (was 200 samples)  ~7 % at 256 Hz
#   gamma band 30-40 Hz (was 30-45)
# Relative power as band / PSD integral (was band / sum of bands) is NOT visible
# here: with nfft = window the default band edges fall on FFT bins and the two
# definitions agree to 1e-12. It is guarded in test-compute_band_power.R, with
# bands that leave a gap.
# Agreement with YASA is about 1e-12 relative, so the tolerance below is 1e-9.

.bp_lcg <- function(n, seed) {       # Park-Miller; 16807 * state < 2^53, exact in doubles
  out <- numeric(n); state <- seed
  for (i in seq_len(n)) {
    state <- (16807 * state) %% 2147483647
    out[i] <- state / 2147483647 - 0.5
  }
  out
}

.bp_signal <- function(fs, n_epochs = 6L) {
  n  <- fs * 30L * n_epochs
  t  <- (seq_len(n) - 1L) / fs
  tp <- 2 * pi
  env <- 1 + 0.6 * sin(tp * t / 1200)
  30 * env * (sin(tp * 1.1 * t + 0.3) + 0.5 * sin(tp * 10.3 * t) +
              0.3 * sin(tp * 2.7 * t + 1)) + 12 * .bp_lcg(n, 12345)
}

.bp_psg <- function(fs, n_epochs = 6L) {
  sig <- .bp_signal(fs, n_epochs)
  spe <- fs * 30L
  structure(
    list(
      edf = list(signals = list(EEG = list(signal = sig))),
      epochs = lapply(seq_len(n_epochs), function(i)
        list(EEG = sig[((i - 1L) * spe + 1L):(i * spe)])),
      n_epochs = n_epochs, epoch_s = 30,
      channel_map = tibble::tibble(label = "EEG", type = "EEG", sample_rate = fs, bad = FALSE)
    ),
    class = "mrpheus_psg"
  )
}

test_that("the synthetic EEG matches the one the band power fixture was made from", {
  # printed by data-raw/generate_yasa_bandpower_parity_fixture.py
  expect_equal(.bp_signal(100)[1:2], c(11.598243405202565, 32.16198145694556),  tolerance = 1e-10)
  expect_equal(.bp_signal(256)[1:2], c(11.598243405202565, 25.274607522481084), tolerance = 1e-10)
})

test_that("compute_band_power equals yasa.bandpower() with default arguments", {
  fix   <- utils::read.csv(testthat::test_path("fixtures", "yasa_bandpower_parity.csv"))
  bands <- c("delta", "theta", "alpha", "sigma", "beta", "gamma")

  for (fs in unique(fix$fs)) {
    f    <- fix[fix$fs == fs, ]
    psg  <- .bp_psg(fs)
    abs_ <- compute_band_power(psg, relative = FALSE)
    rel_ <- compute_band_power(psg, relative = TRUE)
    info <- paste("fs =", fs)

    expect_equal(abs_$epoch, f$epoch, info = info)
    for (b in bands) {
      expect_equal(abs_[[b]], f[[b]],                 tolerance = 1e-9, info = paste(info, b, "absolute"))
      expect_equal(rel_[[b]], f[[paste0("rel_", b)]], tolerance = 1e-9, info = paste(info, b, "relative"))
    }
    expect_equal(abs_$total_power, f$total_power, tolerance = 1e-9, info = paste(info, "total_power"))
  }
})
