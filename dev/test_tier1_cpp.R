# dev/test_tier1_cpp.R
#
# Parity validation for Tier-1 Rcpp functions added in this session:
#   nzc_cpp, petrosian_fd_cpp, hjorth_cpp, stat_features_cpp
#
# Also validates roll_triang_mean_cpp now that the shadowing bug is fixed.
#
# Run after:
#   Rcpp::compileAttributes()
#   devtools::load_all()
#
# Expected output: all checks print [PASS]. Any [FAIL] indicates a parity break.

library(mrpheus)

cat("Tier-1 Rcpp parity validation\n")
cat("==============================\n\n")

# ── Helpers ───────────────────────────────────────────────────────────────────

check <- function(label, cpp_val, r_val, tol = 1e-10) {
  ok <- isTRUE(all.equal(cpp_val, r_val, tolerance = tol, check.names = FALSE))
  cat(sprintf("[%s] %s\n", if (ok) "PASS" else "FAIL", label))
  if (!ok) {
    cat("  C++:", cpp_val, "\n")
    cat("  R:  ", r_val,   "\n")
  }
  invisible(ok)
}

# Pure-R reference implementations (kept for validation only)
.nzc_r        <- function(x) sum(diff(sign(x)) != 0L)
.petrosian_r  <- function(x) {
  N <- length(x); dx <- diff(x)
  nzc_p <- sum((dx[-length(dx)] * dx[-1L]) < 0L)
  log10(N) / (log10(N) + log10(N / (N + 0.4 * nzc_p)))
}
.hjorth_r <- function(x) {
  d1 <- diff(x); d2 <- diff(d1)
  eps <- .Machine$double.eps
  hmob  <- sqrt(var(d1) / (var(x) + eps))
  hcomp <- sqrt(var(d2) / (var(d1) + eps)) / (hmob + eps)
  c(hmob = hmob, hcomp = hcomp)
}
.stat_r <- function(x) {
  mu <- mean(x); sig <- sd(x); eps <- .Machine$double.eps
  if (sig < eps) return(c(std = 0, iqr = 0, skew = 0, kurt = 0))
  zx <- (x - mu) / sig
  c(std = sig, iqr = IQR(x), skew = mean(zx^3), kurt = mean(zx^4) - 3)
}
.roll_triang_r <- function(x, k = 15L) {
  n <- length(x); half <- (k - 1L) %/% 2L
  w_full <- c(seq_len(half + 1L), seq(half, 1L))
  vapply(seq_len(n), function(i) {
    i_start <- max(1L, i - half); i_end <- min(n, i + half)
    w_start <- (i_start - (i - half)) + 1L
    w_sub <- w_full[w_start:(w_start + (i_end - i_start))]
    sum(x[i_start:i_end] * w_sub) / sum(w_sub)
  }, numeric(1))
}

# ── Test signals ──────────────────────────────────────────────────────────────

set.seed(42)
x_eeg  <- rnorm(3000, mean = 5, sd = 50)   # typical 30-s epoch at 100 Hz
x_short <- rnorm(100)                       # short signal
x_flat  <- rep(2.5, 200)                    # constant (tests zero-variance guard)

# Signals with known properties
x_sine  <- sin(seq(0, 10 * pi, length.out = 1000))  # clean sine: predictable nzc
x_step  <- c(rep(-1, 500), rep(1, 500))              # one sign change only

cat("--- nzc_cpp ---\n")
check("random EEG signal",     nzc_cpp(x_eeg),   .nzc_r(x_eeg))
check("short signal",          nzc_cpp(x_short), .nzc_r(x_short))
check("sine wave",             nzc_cpp(x_sine),  .nzc_r(x_sine))
check("step function (nzc=1)", nzc_cpp(x_step),  .nzc_r(x_step))
check("constant (nzc=0)",      nzc_cpp(x_flat),  .nzc_r(x_flat))

cat("\n--- petrosian_fd_cpp ---\n")
check("random EEG signal",   petrosian_fd_cpp(x_eeg),   .petrosian_r(x_eeg))
check("short signal",        petrosian_fd_cpp(x_short), .petrosian_r(x_short))
check("sine wave",           petrosian_fd_cpp(x_sine),  .petrosian_r(x_sine))
check("step function",       petrosian_fd_cpp(x_step),  .petrosian_r(x_step))

cat("\n--- hjorth_cpp ---\n")
h_cpp <- hjorth_cpp(x_eeg); h_r <- .hjorth_r(x_eeg)
check("hmob  (random EEG)", h_cpp["hmob"],  h_r["hmob"])
check("hcomp (random EEG)", h_cpp["hcomp"], h_r["hcomp"])
h2_cpp <- hjorth_cpp(x_sine); h2_r <- .hjorth_r(x_sine)
check("hmob  (sine)",       h2_cpp["hmob"],  h2_r["hmob"])
check("hcomp (sine)",       h2_cpp["hcomp"], h2_r["hcomp"])

cat("\n--- stat_features_cpp ---\n")
s_cpp <- stat_features_cpp(x_eeg); s_r <- .stat_r(x_eeg)
check("std  (random EEG)",  s_cpp["std"],  s_r["std"])
check("iqr  (random EEG)",  s_cpp["iqr"],  s_r["iqr"])
check("skew (random EEG)",  s_cpp["skew"], s_r["skew"])
check("kurt (random EEG)",  s_cpp["kurt"], s_r["kurt"])
s2_cpp <- stat_features_cpp(x_flat); s2_r <- .stat_r(x_flat)
check("constant: all zeros", unname(s2_cpp), unname(s2_r))

cat("\n--- roll_triang_mean_cpp (shadowing bug fix) ---\n")
x_norm <- rnorm(200)
check("roll_triang_mean k=15",
      roll_triang_mean_cpp(x_norm, 15L),
      .roll_triang_r(x_norm, 15L),
      tol = 1e-12)
check("roll_triang_mean k=7",
      roll_triang_mean_cpp(x_norm, 7L),
      .roll_triang_r(x_norm, 7L),
      tol = 1e-12)

cat("\n--- .roll_triang_mean calls Rcpp (not pure-R) ---\n")
# If the shadowing bug were still present, these would be identical in value
# but this confirms only one definition exists
n_defs <- length(grep("\\.roll_triang_mean <- function",
                       deparse(body(mrpheus:::.add_norm_variants)),
                       fixed = TRUE))
# The body check isn't meaningful here; instead confirm no pure-R loop in .roll_triang_mean
fn_body <- deparse(body(mrpheus:::.roll_triang_mean))
uses_rcpp <- any(grepl("roll_triang_mean_cpp", fn_body))
cat(sprintf("[%s] .roll_triang_mean body calls roll_triang_mean_cpp\n",
            if (uses_rcpp) "PASS" else "FAIL"))

cat("\nDone.\n")
