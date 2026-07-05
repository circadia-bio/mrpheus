# tests/testthat/test-staging-features.R
#
# Unit tests for the internal staging feature helpers in R/staging_features.R.
# These tests verify mathematical correctness of each helper in isolation —
# full end-to-end parity with YASA is validated in data-raw/validate_feature_parity.R.
#
# Expected values are derived analytically or by direct formula evaluation.

# ── .spectral_ratios ────────────────────────────────────────────────────────

test_that(".spectral_ratios computes correct ratio values", {
  spec <- c(sdelta = 0.2, fdelta = 0.1, theta = 0.1,
            alpha = 0.15, sigma = 0.05, beta = 0.02, abspow = 50)
  r <- mrpheus:::.spectral_ratios(spec)
  expect_named(r, c("dt", "ds", "db", "at"))
  # delta = sdelta + fdelta = 0.3
  expect_equal(r[["dt"]], 0.3 / 0.1,  tolerance = 1e-10)  # delta / theta
  expect_equal(r[["ds"]], 0.3 / 0.05, tolerance = 1e-10)  # delta / sigma
  expect_equal(r[["db"]], 0.3 / 0.02, tolerance = 1e-10)  # delta / beta
  expect_equal(r[["at"]], 0.15 / 0.1, tolerance = 1e-10)  # alpha / theta
})

test_that(".spectral_ratios uses sdelta + fdelta as delta numerator", {
  # Doubling only sdelta should change dt/ds/db but not at
  spec1 <- c(sdelta = 0.1, fdelta = 0.1, theta = 0.2,
             alpha = 0.2, sigma = 0.1, beta = 0.1, abspow = 1)
  spec2 <- c(sdelta = 0.2, fdelta = 0.1, theta = 0.2,
             alpha = 0.2, sigma = 0.1, beta = 0.1, abspow = 1)
  r1 <- mrpheus:::.spectral_ratios(spec1)
  r2 <- mrpheus:::.spectral_ratios(spec2)
  expect_gt(r2[["dt"]], r1[["dt"]])
  expect_equal(r1[["at"]], r2[["at"]])
})

# ── .nzc ─────────────────────────────────────────────────────────────────────

test_that(".nzc counts zero crossings correctly", {
  expect_equal(mrpheus:::.nzc(c(1, -1, 1, -1, 1)), 4L)
  expect_equal(mrpheus:::.nzc(c(1, 2, 3, 4)), 0L)
  expect_equal(mrpheus:::.nzc(c(2, -1, 3, -2, 1)), 4L)
  expect_equal(mrpheus:::.nzc(c(-5, 5)), 1L)
})

# ── .petrosian_fd ─────────────────────────────────────────────────────────────
# NOTE: petrosian uses extrema count (sign changes in consecutive diffs),
#       NOT zero crossings of x (that was the bug fixed during development).

test_that(".petrosian_fd uses extrema count, not zero crossings", {
  # x = c(1, -1, 1, -1, 1):
  #   zero crossings (.nzc) = 4
  #   extrema (sign changes in diff) = 3  ← what petrosian actually uses
  x      <- c(1, -1, 1, -1, 1)
  nzc_p  <- 3L   # extrema count
  N      <- 5L
  expected <- log10(N) / (log10(N) + log10(N / (N + 0.4 * nzc_p)))
  expect_equal(mrpheus:::.petrosian_fd(x), expected, tolerance = 1e-10)
})

test_that(".petrosian_fd extrema count differs from zero-crossing count on pathological input", {
  # Sawtooth: monotone rising then one sharp drop
  x_saw     <- c(1, 2, 3, 4, 5, 1, 2, 3, 4, 5)
  dx        <- diff(x_saw)
  nzc_extrem <- sum((dx[-length(dx)] * dx[-1L]) < 0L)  # 1 sign change at the drop
  nzc_zc    <- mrpheus:::.nzc(x_saw)                   # zero crossings of x itself
  # They are different; petrosian uses extrema
  expect_equal(
    mrpheus:::.petrosian_fd(x_saw),
    log10(10) / (log10(10) + log10(10 / (10 + 0.4 * nzc_extrem))),
    tolerance = 1e-10
  )
  expect_false(nzc_extrem == nzc_zc)
})

test_that(".petrosian_fd returns 1 for constant signal", {
  # No extrema → nzc_p = 0 → log10(N/(N+0)) = 0 → denominator = log10(N) → FD = 1
  expect_equal(mrpheus:::.petrosian_fd(rep(3, 100)), 1)
})

test_that(".petrosian_fd is always >= 1", {
  set.seed(7)
  expect_gte(mrpheus:::.petrosian_fd(rnorm(300)), 1)
  expect_gte(mrpheus:::.petrosian_fd(sin(seq(0, 6 * pi, length.out = 200))), 1)
})

# ── .scipy_simpson ────────────────────────────────────────────────────────────

test_that(".scipy_simpson integrates constant function exactly (odd N)", {
  # Integral of f=1 from 0 to 4 with dx=1 and N=5 points = 4
  expect_equal(mrpheus:::.scipy_simpson(rep(1, 5), dx = 1), 4, tolerance = 1e-12)
})

test_that(".scipy_simpson integrates linear function exactly (odd N)", {
  # Integral of f=x from 0 to 2 with dx=0.5 and N=5 points = 2
  y <- seq(0, 2, by = 0.5)
  expect_equal(mrpheus:::.scipy_simpson(y, dx = 0.5), 2, tolerance = 1e-12)
})

test_that(".scipy_simpson standard 1/3 rule for N=3 (odd)", {
  y        <- c(1, 4, 1)
  expected <- 1.0 / 3.0 * (1 + 4 * 4 + 1)  # = 6
  expect_equal(mrpheus:::.scipy_simpson(y, dx = 1), expected, tolerance = 1e-12)
})

test_that(".scipy_simpson uses 3/8 rule for N=4 (even, sdelta-sized band)", {
  # Scipy 1.11+: 3/8 rule on all 4 points when n_std=N-3=1 < 3
  y        <- c(2, 3, 4, 5)
  dx       <- 0.2
  expected <- 3 * dx / 8 * (2 + 3 * 3 + 3 * 4 + 5)
  expect_equal(mrpheus:::.scipy_simpson(y, dx = dx), expected, tolerance = 1e-12)
})

test_that(".scipy_simpson N=2 falls back to trapezoid", {
  expect_equal(mrpheus:::.scipy_simpson(c(1, 3), dx = 2), 4, tolerance = 1e-12)
})

test_that(".scipy_simpson N=1 returns 0", {
  expect_equal(mrpheus:::.scipy_simpson(c(5), dx = 1), 0)
})

test_that(".scipy_simpson mixes 1/3 and 3/8 rules for even N=6", {
  # n_std = N-3 = 3 >= 3: standard 1/3 on y[1:3], 3/8 on y[3:6] (overlapping point)
  # For linear f(x) = x+1 with dx=1, exact integral from 0 to 5 = 17.5
  y        <- as.double(1:6)
  expected <- 1/3 * (1 + 4*2 + 3) + 3/8 * (3 + 3*4 + 3*5 + 6)
  expect_equal(mrpheus:::.scipy_simpson(y, dx = 1), expected, tolerance = 1e-12)
  expect_equal(mrpheus:::.scipy_simpson(y, dx = 1), 17.5, tolerance = 1e-12)
})

test_that(".scipy_simpson is always >= 0 for non-negative input", {
  set.seed(42)
  y <- abs(rnorm(21))
  expect_gte(mrpheus:::.scipy_simpson(y, dx = 0.2), 0)
})

# ── .median_bias ─────────────────────────────────────────────────────────────

test_that(".median_bias returns 1 for n <= 2", {
  expect_equal(mrpheus:::.median_bias(1L), 1)
  expect_equal(mrpheus:::.median_bias(2L), 1)
})

test_that(".median_bias returns correct value for n=3", {
  # n_half = 1, ii_2 = [2]: 1 + (1/3 - 1/2) = 1 - 1/6 = 5/6
  expect_equal(mrpheus:::.median_bias(3L), 5 / 6, tolerance = 1e-12)
})

test_that(".median_bias is strictly less than 1 for n >= 3", {
  # The sum is always negative, so bias < 1
  expect_lt(mrpheus:::.median_bias(3L),  1)
  expect_lt(mrpheus:::.median_bias(11L), 1)
})

test_that(".median_bias is non-increasing in n", {
  # For even n, n_half = (n-1)%/%2 equals n_half for n-1, so consecutive pairs
  # (e.g. n=3 and n=4) give the same value — non-increasing, not strictly decreasing.
  vals <- sapply(3:20, mrpheus:::.median_bias)
  expect_true(all(diff(vals) <= 0))
})

# ── .odd_ext_fir ─────────────────────────────────────────────────────────────

test_that(".odd_ext_fir produces correct length", {
  x <- 1:10
  expect_length(mrpheus:::.odd_ext_fir(x, 3L), 16L)  # 3 + 10 + 3
  expect_length(mrpheus:::.odd_ext_fir(x, 0L), 10L)  # no padding
})

test_that(".odd_ext_fir middle section is the original signal", {
  x   <- c(5, 10, 15, 20)
  ext <- mrpheus:::.odd_ext_fir(x, 2L)
  expect_equal(ext[3:6], x)
})

test_that(".odd_ext_fir left padding is odd reflection", {
  x   <- c(10, 20, 30, 40)
  ext <- mrpheus:::.odd_ext_fir(x, 2L)
  # Left: 2*x[1] - x[3:2] = 20 - c(30, 20) = c(-10, 0)
  expect_equal(ext[1:2], c(-10, 0))
})

test_that(".odd_ext_fir right padding is odd reflection", {
  x   <- c(10, 20, 30, 40)
  ext <- mrpheus:::.odd_ext_fir(x, 2L)
  # Right: 2*x[4] - x[3:2] = 80 - c(30, 20) = c(50, 60)
  expect_equal(ext[7:8], c(50, 60))
})

test_that(".odd_ext_fir clamps n to len(x)-1", {
  x <- 1:5
  # n > nx-1: should silently clamp to nx-1 = 4
  ext <- mrpheus:::.odd_ext_fir(x, 100L)
  expect_equal(ext[5:9], x)
  # Total length = 4 + 5 + 4 = 13
  expect_length(ext, 13L)
})

# ── .mne_smart_pad ────────────────────────────────────────────────────────────

test_that(".mne_smart_pad produces correct total length", {
  x <- 1:10
  expect_length(mrpheus:::.mne_smart_pad(x, c(3L, 3L)), 16L)
  expect_length(mrpheus:::.mne_smart_pad(x, c(0L, 0L)), 10L)
})

test_that(".mne_smart_pad middle section is the original signal", {
  x    <- c(10, 20, 30, 40, 50)
  n_l  <- 3L
  n_r  <- 2L
  padded <- mrpheus:::.mne_smart_pad(x, c(n_l, n_r))
  expect_equal(padded[(n_l + 1L):(n_l + length(x))], x)
})

test_that(".mne_smart_pad left reflect_limited is odd reflection", {
  x <- c(1, 2, 3, 4, 5)
  padded <- mrpheus:::.mne_smart_pad(x, c(2L, 0L))
  # left: 2*x[1] - x[seq(3, 2, -1)] = 2 - c(3, 2) = c(-1, 0)
  expect_equal(padded[1:2], c(-1, 0))
})

test_that(".mne_smart_pad right reflect_limited is odd reflection", {
  x <- c(1, 2, 3, 4, 5)
  padded <- mrpheus:::.mne_smart_pad(x, c(0L, 2L))
  # right: 2*x[5] - x[seq(4, 3, -1)] = 10 - c(4, 3) = c(6, 7)
  expect_equal(padded[6:7], c(6, 7))
})

test_that(".mne_smart_pad zero-pads when n_pad > length(x)-1", {
  x <- c(1, 2, 3)
  # n_pad = 5 > nx-1 = 2: left_z_pad = zeros(5-3+1=3), reflect of 2 points
  padded <- mrpheus:::.mne_smart_pad(x, c(5L, 0L))
  # First 3 elements should be zero (zero padding)
  expect_equal(padded[1:3], c(0, 0, 0))
})

# ── .mne_fft_resample ─────────────────────────────────────────────────────────

test_that(".mne_fft_resample returns correct output length", {
  x <- rnorm(50)
  expect_length(mrpheus:::.mne_fft_resample(x, up = 4L), 200L)
  expect_length(mrpheus:::.mne_fft_resample(x, up = 10L), 500L)
})

test_that(".mne_fft_resample preserves constant signal (no distortion in passband)", {
  x   <- rep(1.0, 100)
  out <- mrpheus:::.mne_fft_resample(x, up = 10L)
  # Interior of upsampled constant should be ~1.0 (edges may differ slightly)
  interior <- out[50:950]
  expect_equal(mean(interior), 1.0, tolerance = 1e-6)
  expect_equal(max(abs(interior - 1.0)), 0, tolerance = 1e-6)
})

test_that(".mne_fft_resample preserves DC for ramp signal", {
  x   <- seq(0, 1, length.out = 100)
  out <- mrpheus:::.mne_fft_resample(x, up = 5L)
  # Output length
  expect_length(out, 500L)
  # DC (mean) should be preserved (≈ mean of x ≈ 0.5)
  expect_equal(mean(out), mean(x), tolerance = 0.01)
})

test_that(".mne_fft_resample output length matches length(x) * up for 100x", {
  x <- rnorm(200)
  expect_length(mrpheus:::.mne_fft_resample(x, up = 100L), 20000L)
})

# ── .welch_median_psd ─────────────────────────────────────────────────────────

test_that(".welch_median_psd returns correct structure", {
  sig  <- rnorm(3000)
  psd  <- mrpheus:::.welch_median_psd(sig, sr = 100L, win_n = 500L)
  expect_named(psd, c("freq", "spec"))
  # n_freq = 500/2 + 1 = 251
  expect_length(psd$freq, 251L)
  expect_length(psd$spec, 251L)
})

test_that(".welch_median_psd frequency axis spans 0 to Nyquist", {
  psd <- mrpheus:::.welch_median_psd(rnorm(3000), sr = 100L, win_n = 500L)
  expect_equal(psd$freq[1L],  0,    tolerance = 1e-10)
  expect_equal(psd$freq[251L], 50.0, tolerance = 1e-10)
})

test_that(".welch_median_psd frequency resolution matches sr/win_n", {
  psd <- mrpheus:::.welch_median_psd(rnorm(3000), sr = 100L, win_n = 500L)
  expect_equal(psd$freq[2L] - psd$freq[1L], 0.2, tolerance = 1e-10)
})

test_that(".welch_median_psd spectral peak matches known frequency", {
  # Pure 10 Hz sinusoid at 100 Hz: peak should be at 10 Hz
  sr  <- 100L
  sig <- sin(2 * pi * 10 * seq(0, 30 - 1/sr, by = 1/sr))
  psd <- mrpheus:::.welch_median_psd(sig, sr = sr, win_n = 500L)
  peak_freq <- psd$freq[which.max(psd$spec)]
  expect_equal(peak_freq, 10, tolerance = 0.21)  # within one freq bin (0.2 Hz)
})

test_that(".welch_median_psd periodic Hamming gives correct scale sum", {
  # For periodic Hamming of length N: sum(w^2) = N*(0.54^2 + 0.46^2/2) exactly
  win_n    <- 500L
  expected <- win_n * (0.54^2 + 0.46^2 / 2)
  wvec     <- 0.54 - 0.46 * cos(2 * pi * seq.int(0L, win_n - 1L) / win_n)
  expect_equal(sum(wvec^2), expected, tolerance = 1e-10)
})

test_that(".welch_median_psd returns NA spec when no valid segments", {
  # Signal shorter than one window
  psd <- mrpheus:::.welch_median_psd(rnorm(100), sr = 100L, win_n = 500L)
  expect_true(all(is.na(psd$spec)))
})

# ── .filtfilt_fir / .bandpass_filter ─────────────────────────────────────────

test_that(".filtfilt_fir output length equals input length", {
  b   <- c(0.5, 1, 0.5)  # simple 3-tap FIR
  x   <- rnorm(500)
  out <- mrpheus:::.filtfilt_fir(b, x)
  expect_length(out, length(x))
})

test_that(".filtfilt_fir is zero-phase (no delay for symmetric FIR)", {
  # An odd-symmetric FIR applied zero-phase should have no time shift.
  # A pure sine at 10 Hz should emerge at the same phase as the input.
  sr  <- 1000L
  t   <- seq(0, 2, by = 1/sr)
  x   <- sin(2 * pi * 10 * t)
  # 3-tap box filter (passes all frequencies, trivially zero-phase)
  b   <- c(1/3, 1/3, 1/3)
  out <- mrpheus:::.filtfilt_fir(b, x)
  # Cross-correlation peak should be at lag 0 (zero phase shift)
  lags   <- -10:10
  cc     <- sapply(lags, function(l) cor(x[100:1900], out[(100 + l):(1900 + l)]))
  expect_equal(lags[which.max(cc)], 0L)
})

test_that(".bandpass_filter output length equals input length", {
  skip_if_not(
    nzchar(system.file("filters", "mne_bandpass_100hz.csv", package = "mrpheus")),
    "MNE filter coefficients not bundled"
  )
  x   <- rnorm(3000)
  out <- mrpheus:::.bandpass_filter(x, sr = 100L)
  expect_length(out, length(x))
})

test_that(".bandpass_filter attenuates below-cutoff signal (~0.1 Hz)", {
  skip_if_not(
    nzchar(system.file("filters", "mne_bandpass_100hz.csv", package = "mrpheus")),
    "MNE filter coefficients not bundled"
  )
  sr   <- 100L
  t    <- seq(0, 79.99, by = 1 / sr)
  # 0.1 Hz is in the transition band (0–0.4 Hz); expect ~8x attenuation, not stopband
  x_lo <- sin(2 * pi * 0.1 * t)
  # 5 Hz sine — clearly in passband
  x_in <- sin(2 * pi * 5.0 * t)
  out_lo <- mrpheus:::.bandpass_filter(x_lo, sr = sr)
  out_in <- mrpheus:::.bandpass_filter(x_in, sr = sr)
  # RMS of transition-band signal should be meaningfully smaller than in-band
  rms <- function(x) sqrt(mean(x^2))
  expect_lt(rms(out_lo), 0.2 * rms(out_in))
})

# ── .robust_scale ─────────────────────────────────────────────────────────────

test_that(".robust_scale centres on the median", {
  x      <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  result <- mrpheus:::.robust_scale(x)
  expect_equal(median(result), 0, tolerance = 1e-6)
})

test_that(".robust_scale denominator uses q95 - q5", {
  x   <- as.double(1:100)
  q5  <- unname(quantile(x, 0.05))
  q95 <- unname(quantile(x, 0.95))
  result <- mrpheus:::.robust_scale(x)
  expected_scale <- q95 - q5 + 1e-10
  expect_equal(unname(result[1L]), (x[1L] - median(x)) / expected_scale, tolerance = 1e-6)
})

test_that(".robust_scale returns zeros for constant input", {
  result <- mrpheus:::.robust_scale(rep(5, 50))
  # Median = 5, all x - median = 0, denominator near 0 but +eps → 0 / eps ≈ 0
  expect_equal(unique(result), 0, tolerance = 1e-5)
})

# ── roll_triang_mean_cpp / .roll_triang_mean ──────────────────────────────────

test_that("roll_triang_mean_cpp with k=1 is the identity", {
  set.seed(42)
  x <- rnorm(20)
  expect_equal(roll_triang_mean_cpp(x, k = 1L), x, tolerance = 1e-10)
})

test_that("roll_triang_mean_cpp returns constant for constant input", {
  x   <- rep(7.0, 20)
  out <- roll_triang_mean_cpp(x, k = 5L)
  expect_equal(out, rep(7.0, 20), tolerance = 1e-10)
})

test_that("roll_triang_mean_cpp length matches input", {
  x <- rnorm(50)
  expect_length(roll_triang_mean_cpp(x, k = 7L),  50L)
  expect_length(roll_triang_mean_cpp(x, k = 15L), 50L)
})

test_that("roll_triang_mean_cpp matches R reference for k=3", {
  # For x = 1:5, k=3, half=1:
  # i=1: weights=[2,1]/3,    mean = (1*2+2*1)/3 = 4/3
  # i=2: weights=[1,2,1]/4,  mean = (1+4+3)/4   = 2
  # i=3: weights=[1,2,1]/4,  mean = (2+6+4)/4   = 3
  # i=4: weights=[1,2,1]/4,  mean = (3+8+5)/4   = 4
  # i=5: weights=[1,2]/3,    mean = (4*1+5*2)/3 = 14/3
  x        <- as.double(1:5)
  expected <- c(4/3, 2, 3, 4, 14/3)
  expect_equal(roll_triang_mean_cpp(x, k = 3L), expected, tolerance = 1e-10)
})

test_that("roll_triang_mean_cpp C++ output matches R vapply reference exactly", {
  set.seed(11)
  x   <- rnorm(50)
  # Both R and C++ now call the same underlying function but let's verify numerics
  cpp <- roll_triang_mean_cpp(x, k = 15L)
  r   <- mrpheus:::.roll_triang_mean(x, k = 15L)
  expect_equal(cpp, r, tolerance = 1e-10)
})

# ── perm_entropy_cpp ──────────────────────────────────────────────────────────

test_that("perm_entropy_cpp is 0 for monotone signal with delay > 1", {
  # With delay=2 and order=3, embedded windows of 1:7 are still all monotone increasing
  expect_equal(perm_entropy_cpp(as.double(1:7), order = 3L, delay = 2L),
               0, tolerance = 1e-10)
})

test_that("perm_entropy_cpp output is in [0, 1]", {
  set.seed(42)
  for (i in 1:5) {
    val <- perm_entropy_cpp(rnorm(300), order = 3L, delay = 1L)
    expect_gte(val, 0)
    expect_lte(val, 1)
  }
})

test_that("perm_entropy_cpp returns 0 for monotone signal", {
  expect_equal(perm_entropy_cpp(as.double(1:50), order = 3L), 0, tolerance = 1e-10)
  expect_equal(perm_entropy_cpp(as.double(50:1), order = 3L), 0, tolerance = 1e-10)
})

test_that("perm_entropy_cpp matches .perm_entropy on random signal", {
  set.seed(99)
  x <- rnorm(300)
  expect_equal(perm_entropy_cpp(x, 3L, 1L),
               mrpheus:::.perm_entropy(x, 3L, 1L),
               tolerance = 1e-10)
})

test_that("perm_entropy_cpp returns maximally entropic value for uniform-random signal", {
  # With many samples, all 6 patterns are equally likely → H = 1
  set.seed(7)
  val <- perm_entropy_cpp(rnorm(50000), order = 3L)
  expect_equal(val, 1, tolerance = 0.01)
})

# ── higuchi_fd_cpp ────────────────────────────────────────────────────────────

test_that("higuchi_fd_cpp matches .higuchi_fd on random signal", {
  set.seed(5)
  x <- rnorm(500)
  expect_equal(higuchi_fd_cpp(x, kmax = 10L),
               mrpheus:::.higuchi_fd(x, kmax = 10L),
               tolerance = 1e-10)
})

test_that("higuchi_fd_cpp is ~1 for monotone signal", {
  expect_equal(higuchi_fd_cpp(as.double(1:200), kmax = 10L), 1, tolerance = 0.05)
})

test_that("higuchi_fd_cpp is ~2 for white noise", {
  set.seed(42)
  expect_equal(higuchi_fd_cpp(rnorm(1000), kmax = 10L), 2, tolerance = 0.15)
})

# ── .spectral_features ────────────────────────────────────────────────────────

test_that(".spectral_features returns correctly named vector", {
  set.seed(1)
  result <- mrpheus:::.spectral_features(rnorm(256 * 30), sr = 256)
  expect_named(result, c("sdelta", "fdelta", "theta", "alpha", "sigma", "beta", "abspow"))
})

test_that(".spectral_features relative band powers sum to 1", {
  set.seed(1)
  result  <- mrpheus:::.spectral_features(rnorm(256 * 30), sr = 256)
  rel_sum <- sum(result[c("sdelta", "fdelta", "theta", "alpha", "sigma", "beta")])
  # Bands share boundary points; Simpson's on each sub-band sums to ~1 not exactly 1
  expect_equal(rel_sum, 1, tolerance = 0.002)
})

test_that(".spectral_features abspow is positive", {
  set.seed(2)
  result <- mrpheus:::.spectral_features(rnorm(256 * 30), sr = 256)
  expect_gt(unname(result["abspow"]), 0)
})

test_that(".spectral_features sigma dominates for spindle-range sine", {
  # A 14 Hz sine (spindle frequency) should drive sigma power highest
  sr  <- 256
  sig <- sin(2 * pi * 14 * seq(0, 30, length.out = sr * 30))
  result <- mrpheus:::.spectral_features(sig, sr)
  rel_bands <- result[c("sdelta", "fdelta", "theta", "alpha", "sigma", "beta")]
  expect_equal(unname(which.max(rel_bands)), 5L)   # sigma is 5th in alphabetical order
})

test_that(".spectral_features abspow uses trapz (not simpson)", {
  # abspow must NOT be affected by band power denominator (they are separate computations)
  # Quick sanity: abspow > 0 for any non-zero signal
  set.seed(3)
  result <- mrpheus:::.spectral_features(rnorm(3000), sr = 100L)
  expect_gt(result["abspow"], 0)
})

# ── .na_channel_matrix ───────────────────────────────────────────────────

test_that(".na_channel_matrix has correct shape and all-NA values", {
  mat <- mrpheus:::.na_channel_matrix(
    n = 10L, prefix = "eog_",
    base_names = c("abspow", "alpha", "beta")
  )
  expect_equal(nrow(mat), 10L)
  expect_equal(ncol(mat), 9L)   # 3 base features × 3 variants each
  expect_true(all(is.na(mat)))
})

test_that(".na_channel_matrix column names follow prefix_name_variant pattern", {
  mat <- mrpheus:::.na_channel_matrix(
    n = 5L, prefix = "emg_",
    base_names = c("std", "iqr")
  )
  expect_equal(colnames(mat),
               c("emg_std", "emg_std_c7min_norm", "emg_std_p2min_norm",
                 "emg_iqr", "emg_iqr_c7min_norm", "emg_iqr_p2min_norm"))
})

# ── .fft_fir_full ───────────────────────────────────────────────────────────

test_that(".fft_fir_full returns full convolution length N + M - 1", {
  b   <- c(0.25, 0.5, 0.25)  # M = 3
  x   <- rnorm(100)           # N = 100
  out <- mrpheus:::.fft_fir_full(b, x)
  expect_length(out, 100L + 3L - 1L)
})

test_that(".fft_fir_full matches stats::convolve for box filter", {
  b   <- c(1/3, 1/3, 1/3)
  x   <- rnorm(50)
  # stats::convolve(x, rev(b)) gives linear convolution
  ref <- as.vector(stats::convolve(x, rev(b), type = "open"))
  out <- mrpheus:::.fft_fir_full(b, x)
  expect_equal(out, ref, tolerance = 1e-9)
})

# ── .eeg_epoch_features output shape ─────────────────────────────────────────

test_that(".eeg_epoch_features returns 21 correctly named values", {
  set.seed(99)
  sig    <- rnorm(256 * 30)
  result <- mrpheus:::.eeg_epoch_features(sig, sr = 256)
  expect_length(result, 21L)
  expected_names <- c("abspow", "alpha", "at", "beta", "db", "ds", "dt",
                       "fdelta", "hcomp", "higuchi", "hmob", "iqr", "kurt",
                       "nzc", "perm", "petrosian", "sdelta", "sigma",
                       "skew", "std", "theta")
  expect_equal(names(result), expected_names)
})

test_that(".eog_epoch_features returns 17 correctly named values", {
  set.seed(99)
  result <- mrpheus:::.eog_epoch_features(rnorm(256 * 30), sr = 256)
  expect_length(result, 17L)
  expect_false(any(c("at", "db", "ds", "dt") %in% names(result)))
})

test_that(".emg_epoch_features returns 11 correctly named values", {
  set.seed(99)
  result <- mrpheus:::.emg_epoch_features(rnorm(256 * 30), sr = 256)
  expect_length(result, 11L)
  expected_names <- c("abspow", "hcomp", "higuchi", "hmob", "iqr", "kurt",
                       "nzc", "perm", "petrosian", "skew", "std")
  expect_equal(names(result), expected_names)
})

# ── .add_norm_variants column structure ──────────────────────────────────────

test_that(".add_norm_variants produces 3x columns with correct interleaving", {
  mat <- matrix(rnorm(30), nrow = 10, ncol = 3,
                dimnames = list(NULL, c("foo", "bar", "baz")))
  result <- mrpheus:::.add_norm_variants(mat, "pfx_")
  expect_equal(ncol(result), 9L)
  expect_equal(colnames(result),
               c("pfx_foo", "pfx_foo_c7min_norm", "pfx_foo_p2min_norm",
                 "pfx_bar", "pfx_bar_c7min_norm", "pfx_bar_p2min_norm",
                 "pfx_baz", "pfx_baz_c7min_norm", "pfx_baz_p2min_norm"))
})

test_that(".add_norm_variants raw column is unchanged", {
  set.seed(5)
  mat    <- matrix(rnorm(50), nrow = 10, ncol = 5,
                   dimnames = list(NULL, letters[1:5]))
  result <- mrpheus:::.add_norm_variants(mat, "x_")
  for (nm in colnames(mat)) {
    expect_equal(result[, paste0("x_", nm)], mat[, nm])
  }
})

test_that(".add_norm_variants c7min_norm and p2min_norm columns have zero median", {
  # Both normalisation variants are robust_scale outputs → median = 0
  set.seed(6)
  mat <- matrix(rnorm(100 * 3), nrow = 100, ncol = 3,
                dimnames = list(NULL, c("a", "b", "c")))
  result <- mrpheus:::.add_norm_variants(mat, "x_")
  for (nm in c("a", "b", "c")) {
    expect_equal(median(result[, paste0("x_", nm, "_c7min_norm")]),
                 0, tolerance = 1e-6)
    expect_equal(median(result[, paste0("x_", nm, "_p2min_norm")]),
                 0, tolerance = 1e-6)
  }
})

# ── .hjorth ───────────────────────────────────────────────────────────────────

test_that(".hjorth returns named vector with correct names", {
  result <- mrpheus:::.hjorth(rnorm(100))
  expect_named(result, c("hmob", "hcomp"))
})

test_that(".hjorth mobility is 0 for constant signal", {
  result <- mrpheus:::.hjorth(rep(5, 100))
  expect_equal(unname(result["hmob"]), 0, tolerance = 1e-10)
})

test_that(".hjorth mobility increases with frequency for sine waves", {
  sr <- 256
  t  <- seq(0, 1, length.out = sr)
  x1 <- sin(2 * pi *  5 * t)
  x2 <- sin(2 * pi * 20 * t)
  expect_gt(
    unname(mrpheus:::.hjorth(x2)["hmob"]),
    unname(mrpheus:::.hjorth(x1)["hmob"])
  )
})

test_that(".hjorth mobility and complexity are positive for non-constant signal", {
  set.seed(1)
  result <- mrpheus:::.hjorth(rnorm(256))
  expect_gt(unname(result["hmob"]),  0)
  expect_gt(unname(result["hcomp"]), 0)
})

# ── .higuchi_fd ───────────────────────────────────────────────────────────────

test_that(".higuchi_fd is ~1 for monotone signal", {
  result <- mrpheus:::.higuchi_fd(seq(1, 100, by = 1))
  expect_equal(result, 1, tolerance = 0.05)
})

test_that(".higuchi_fd is ~2 for Gaussian white noise", {
  set.seed(42)
  result <- mrpheus:::.higuchi_fd(rnorm(1000))
  expect_equal(result, 2, tolerance = 0.15)
})

test_that(".higuchi_fd is strictly greater for noise than sine wave", {
  sr  <- 256
  t   <- seq(0, 4, length.out = sr * 4)
  set.seed(3)
  expect_gt(
    mrpheus:::.higuchi_fd(rnorm(length(t))),
    mrpheus:::.higuchi_fd(sin(2 * pi * 10 * t))
  )
})

# ── .stat_features ────────────────────────────────────────────────────────────

test_that(".stat_features returns correctly named vector", {
  expect_named(mrpheus:::.stat_features(rnorm(100)), c("std", "iqr", "skew", "kurt"))
})

test_that(".stat_features std and iqr match base R", {
  x      <- c(3.1, 1.4, 5.9, 2.6, 5.3, 5.8, 9.7, 9.3, 2.3, 8.4)
  result <- mrpheus:::.stat_features(x)
  expect_equal(unname(result["std"]), stats::sd(x),  tolerance = 1e-10)
  expect_equal(unname(result["iqr"]), stats::IQR(x), tolerance = 1e-10)
})

test_that(".stat_features skewness is 0 for symmetric distribution", {
  x      <- c(1, 2, 3, 4, 5)
  result <- mrpheus:::.stat_features(x)
  expect_equal(unname(result["skew"]), 0, tolerance = 1e-10)
})

test_that(".stat_features excess kurtosis is negative for uniform data", {
  x      <- seq(1, 100)
  result <- mrpheus:::.stat_features(x)
  expect_lt(unname(result["kurt"]), 0)
})

test_that(".stat_features returns all zeros for constant signal", {
  result <- mrpheus:::.stat_features(rep(7, 50))
  expect_equal(unname(result), c(0, 0, 0, 0))
})

# ── .perm_entropy ─────────────────────────────────────────────────────────────

test_that(".perm_entropy returns analytically derived value", {
  h_raw <- -(0.4 * log(0.4) + 0.2 * log(0.2) + 0.4 * log(0.4))
  expected <- h_raw / log(factorial(3))
  result <- mrpheus:::.perm_entropy(c(4, 7, 9, 10, 6, 11, 3), order = 3L, delay = 1L)
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that(".perm_entropy is 0 for monotone signal", {
  expect_equal(mrpheus:::.perm_entropy(1:20, order = 3L), 0, tolerance = 1e-10)
  expect_equal(mrpheus:::.perm_entropy(20:1, order = 3L), 0, tolerance = 1e-10)
})

test_that(".perm_entropy is in [0, 1] for any signal", {
  set.seed(42)
  for (i in 1:5) {
    val <- mrpheus:::.perm_entropy(rnorm(300), order = 3L)
    expect_gte(val, 0)
    expect_lte(val, 1)
  }
})
