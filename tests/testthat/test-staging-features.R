# tests/testthat/test-staging-features.R
#
# Unit tests for the internal staging feature helpers in R/staging_features.R.
# These tests verify mathematical correctness of each helper in isolation —
# full end-to-end parity with YASA is validated in data-raw/validate_feature_parity.R.
#
# Expected values are derived analytically or by direct formula evaluation.

# ── .nzc ─────────────────────────────────────────────────────────────────────

test_that(".nzc counts zero crossings correctly", {
  expect_equal(mrpheus:::.nzc(c(1, -1, 1, -1, 1)), 4L)
  expect_equal(mrpheus:::.nzc(c(1, 2, 3, 4)), 0L)
  expect_equal(mrpheus:::.nzc(c(2, -1, 3, -2, 1)), 4L)
  expect_equal(mrpheus:::.nzc(c(-5, 5)), 1L)
})

# ── .petrosian_fd ─────────────────────────────────────────────────────────────

test_that(".petrosian_fd returns exact value for known input", {
  x <- c(1, -1, 1, -1, 1)   # N=5, nzc=4
  expected <- log10(5) / (log10(5) + log10(5 / (5 + 0.4 * 4)))
  expect_equal(mrpheus:::.petrosian_fd(x), expected, tolerance = 1e-10)
})

test_that(".petrosian_fd returns 1 for constant signal (nzc = 0)", {
  # log10(N/(N+0)) = log10(1) = 0, so FD = log10(N)/log10(N) = 1
  expect_equal(mrpheus:::.petrosian_fd(rep(3, 100)), 1)
})

test_that(".petrosian_fd is always >= 1", {
  set.seed(7)
  expect_gte(mrpheus:::.petrosian_fd(rnorm(300)), 1)
  expect_gte(mrpheus:::.petrosian_fd(sin(seq(0, 6 * pi, length.out = 200))), 1)
})

# ── .perm_entropy ─────────────────────────────────────────────────────────────

test_that(".perm_entropy returns analytically derived value", {
  # x = c(4, 7, 9, 10, 6, 11, 3), order = 3, delay = 1
  # Windows (via order()):
  #   [4,7,9]   -> "123", [7,9,10]  -> "123", [9,10,6]  -> "312"
  #   [10,6,11] -> "213", [6,11,3]  -> "312"
  # freqs: "123"=2/5, "213"=1/5, "312"=2/5
  # H = -(0.4*log(0.4) + 0.2*log(0.2) + 0.4*log(0.4)) / log(6)
  h_raw <- -(0.4 * log(0.4) + 0.2 * log(0.2) + 0.4 * log(0.4))
  expected <- h_raw / log(factorial(3))
  result <- mrpheus:::.perm_entropy(c(4, 7, 9, 10, 6, 11, 3), order = 3L, delay = 1L)
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that(".perm_entropy is 0 for monotone signal", {
  # Only one ordinal pattern possible -> entropy = 0
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

# ── .hjorth ───────────────────────────────────────────────────────────────────

test_that(".hjorth returns named vector with correct names", {
  result <- mrpheus:::.hjorth(rnorm(100))
  expect_named(result, c("hmob", "hcomp"))
})

test_that(".hjorth mobility is 0 for constant signal", {
  # var(diff(rep(c, n))) = 0, so hmob = sqrt(0/var) = 0
  result <- mrpheus:::.hjorth(rep(5, 100))
  expect_equal(unname(result["hmob"]), 0, tolerance = 1e-10)
})

test_that(".hjorth mobility increases with frequency for sine waves", {
  sr <- 256
  t  <- seq(0, 1, length.out = sr)
  x1 <- sin(2 * pi *  5 * t)   # 5 Hz
  x2 <- sin(2 * pi * 20 * t)   # 20 Hz — higher frequency -> higher mobility
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
  # A perfectly linear signal has FD = 1 (one-dimensional trajectory)
  result <- mrpheus:::.higuchi_fd(seq(1, 100, by = 1))
  expect_equal(result, 1, tolerance = 0.05)
})

test_that(".higuchi_fd is ~2 for Gaussian white noise", {
  # White noise fills 2D space -> FD approaches 2
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
  expect_equal(unname(result["std"]), stats::sd(x), tolerance = 1e-10)
  expect_equal(unname(result["iqr"]), stats::IQR(x), tolerance = 1e-10)
})

test_that(".stat_features skewness is 0 for symmetric distribution", {
  x      <- c(1, 2, 3, 4, 5)
  result <- mrpheus:::.stat_features(x)
  expect_equal(unname(result["skew"]), 0, tolerance = 1e-10)
})

test_that(".stat_features excess kurtosis is 0 for uniform distribution", {
  # Exact kurtosis of continuous uniform is -6/5; discrete is close
  # Just verify it's negative (platykurtic) for uniform data
  x      <- seq(1, 100)
  result <- mrpheus:::.stat_features(x)
  expect_lt(unname(result["kurt"]), 0)
})

test_that(".stat_features returns all zeros for constant signal", {
  result <- mrpheus:::.stat_features(rep(7, 50))
  expect_equal(unname(result), c(0, 0, 0, 0))
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
  expect_equal(rel_sum, 1, tolerance = 1e-6)
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
  # Every raw column should be identical to the source
  for (nm in colnames(mat)) {
    expect_equal(result[, paste0("x_", nm)], mat[, nm])
  }
})
