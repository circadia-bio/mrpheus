# R/staging_features.R
#
# Internal feature extraction for AASM sleep staging.
# All functions are package-private (prefixed with `.`) — not exported.
#
# Reproduces the YASA feature pipeline in pure R using base R, gsignal,
# pracma, and zoo. The 149-feature matrix must match the Python implementation
# exactly; parity is validated in data-raw/validate_feature_parity.R.
#
# Reference: Vallat & Walker (2021), eLife 10:e70092.

# ── Band definitions ─────────────────────────────────────────────────────────
# NOTE: band definitions confirmed against YASA SleepStaging.fit() source.

.FREQ_BROAD <- c(0.4, 30)

.STAGING_BANDS <- list(
  sdelta = c(0.4,  1),
  fdelta = c(1,    4),
  theta  = c(4,    8),
  alpha  = c(8,   12),
  sigma  = c(12,  16),
  beta   = c(16,  30)
)

# ── Pre-filter helpers ────────────────────────────────────────────────────────
# YASA resamples to 100 Hz then calls MNE filter_data(l_freq=0.4, h_freq=30).
# The exact 825-tap Hamming FIR coefficients are bundled in inst/filters/
# (generated once by data-raw/extract_mne_filter.py, no Python at runtime).
#
# Applied here with a zero-phase FFT filtfilt matching
# scipy.signal.filtfilt(b, [1], x, padtype='odd').
# Adapted from dynR::bandpass_filter (same scipy-parity approach).

# Odd-reflection padding — matches scipy padtype='odd'.
.odd_ext_fir <- function(x, n) {
  nx <- length(x)
  n  <- min(n, nx - 1L)
  c(
    2 * x[1L] - x[seq.int(n + 1L, 2L,      by = -1L)],
    x,
    2 * x[nx]  - x[seq.int(nx - 1L, nx - n, by = -1L)]
  )
}

# Single causal FIR pass via FFT, returns FULL convolution (length N+M-1).
.fft_fir_full <- function(b, x) {
  M     <- length(b)
  N     <- length(x)
  n_fft <- 2L ^ ceiling(log2(N + M - 1L))
  B     <- stats::fft(c(b, numeric(n_fft - M)))
  X     <- stats::fft(c(x, numeric(n_fft - N)))
  Re(stats::fft(X * B, inverse = TRUE))[seq_len(N + M - 1L)] / n_fft
}

# Zero-phase FIR: single FFT pass then remove the linear group delay (M-1)/2.
# Matches MNE's filter_data which applies a symmetric FIR as a single
# zero-phase pass, NOT as filtfilt (which would square the response).
.filtfilt_fir <- function(b, x) {
  M     <- length(b)
  delay <- (M - 1L) %/% 2L      # group delay of symmetric FIR
  pad   <- 3L * delay
  ext   <- .odd_ext_fir(x, pad)
  # Full convolution, then shift by delay to remove group delay
  y_full     <- .fft_fir_full(b, ext)
  y_nodelay  <- y_full[seq.int(delay + 1L, delay + length(ext))]
  # Strip padding to recover the filtered signal
  y_nodelay[seq.int(pad + 1L, pad + length(x))]
}

.bandpass_filter <- function(sig, sr) {
  coef_path <- system.file("filters", "mne_bandpass_100hz.csv",
                            package = "mrpheus")
  if (!nzchar(coef_path))
    cli::cli_abort("Filter coefficients not found. Run data-raw/extract_mne_filter.py.")
  b <- scan(coef_path, quiet = TRUE)
  .filtfilt_fir(b, sig)
}

# ── Spectral helpers ──────────────────────────────────────────────────────────

# Composite Simpson's rule matching scipy.integrate.simpson.
# - Odd  N: standard 1/3 rule.
# - Even N: scipy 1.11+ approach — standard 1/3 on first N-3 points,
#            Simpson's 3/8 rule on the last 4 points.
.scipy_simpson <- function(y, dx) {
  N <- length(y)
  if (N < 2L) return(0)
  if (N == 2L) return(dx * (y[1L] + y[2L]) / 2)   # degenerate: trapezoid

  if (N %% 2L == 1L) {
    # Odd N: standard composite 1/3 rule
    # coefficients: 1 4 2 4 2 ... 4 1
    coef      <- rep(2, N)
    coef[seq(2L, N - 1L, by = 2L)] <- 4
    coef[c(1L, N)] <- 1
    return(dx / 3 * sum(coef * y))
  } else {
    # Even N: standard 1/3 on first N-3 points (odd count) +
    #          3/8 rule on last 4 points
    n_std   <- N - 3L
    s_total <- 0

    if (n_std >= 3L) {
      coef <- rep(2, n_std)
      coef[seq(2L, n_std - 1L, by = 2L)] <- 4
      coef[c(1L, n_std)] <- 1
      s_total <- dx / 3 * sum(coef * y[seq_len(n_std)])
    }

    # Simpson's 3/8 on last 4 points
    s_total + 3 * dx / 8 * (y[N - 3L] + 3 * y[N - 2L] + 3 * y[N - 1L] + y[N])
  }
}

# Bias correction factor for the median of a chi-squared distribution.
# scipy.signal.welch(average='median') divides by this before returning.
# Matches scipy.signal.spectral._median_bias(n).
.median_bias <- function(n) {
  n_half <- (n - 1L) %/% 2L
  if (n_half < 1L) return(1)
  ii_2 <- 2 * seq_len(n_half)
  1 + sum(1 / (ii_2 + 1) - 1 / ii_2)
}

# Compute Welch PSD with bias-corrected median averaging over segments.
# Matches YASA/scipy.signal.welch(average='median').
.welch_median_psd <- function(sig, sr, win_n, overlap = 0.5) {
  wvec   <- gsignal::hamming(win_n)
  hop    <- max(1L, as.integer(win_n * (1 - overlap)))
  n_freq <- win_n %/% 2L + 1L
  freqs  <- seq(0, sr / 2, length.out = n_freq)
  starts <- seq(1L, length(sig) - win_n + 1L, by = hop)

  if (length(starts) == 0L)
    return(list(freq = freqs, spec = rep(NA_real_, n_freq)))

  scale <- sum(wvec ^ 2) * sr

  pgrams <- vapply(starts, function(s) {
    ft  <- stats::fft(sig[s:(s + win_n - 1L)] * wvec)[seq_len(n_freq)]
    psd <- Mod(ft) ^ 2 / scale
    psd[seq(2L, n_freq - 1L)] <- 2 * psd[seq(2L, n_freq - 1L)]
    psd
  }, numeric(n_freq))

  bias    <- .median_bias(length(starts))
  psd_med <- apply(pgrams, 1L, stats::median) / bias
  list(freq = freqs, spec = psd_med)
}

# Compute Welch PSD and return named vector of band features:
# sdelta, fdelta, theta, alpha, sigma, beta (relative) + abspow.
# Matches YASA: 5-second Hamming window, median averaging.
# Band powers: relative, computed via bandpower_from_psd_ndarray (simpson, inclusive bounds).
# abspow: trapezoid over freq_broad (matches YASA's separate trapezoid call).
.spectral_features <- function(sig, sr) {
  win <- as.integer(5L * sr)
  psd <- .welch_median_psd(sig, sr, win_n = win)
  dx  <- psd$freq[2L] - psd$freq[1L]   # frequency resolution

  # Band powers via Simpson's integration with INCLUSIVE upper bound.
  # Matches YASA's bandpower_from_psd_ndarray (freqs >= b0 & freqs <= b1).
  bp <- vapply(.STAGING_BANDS, function(b) {
    idx <- psd$freq >= b[1] & psd$freq <= b[2]   # inclusive
    if (!any(idx)) return(NA_real_)
    .scipy_simpson(psd$spec[idx], dx)
  }, numeric(1))

  # Total power: Simpson's over full freq_broad range (denominator for relative powers).
  # Matches YASA's: total_power = simpson(psd_trimmed, dx=res)
  idx_broad_full <- psd$freq >= .FREQ_BROAD[1] & psd$freq <= .FREQ_BROAD[2]
  total_power    <- .scipy_simpson(psd$spec[idx_broad_full], dx)
  rel            <- bp / total_power

  # abspow: trapezoid over freq_broad, matching YASA's separate trapezoid call
  abspow <- pracma::trapz(psd$freq[idx_broad_full], psd$spec[idx_broad_full])

  c(rel, abspow = abspow)
}

# Spectral ratio features (EEG only).
# Numerator uses full delta (sdelta + fdelta = 0.4–4 Hz), matching YASA.
.spectral_ratios <- function(spec) {
  delta <- spec[["sdelta"]] + spec[["fdelta"]]
  c(
    dt = delta           / spec[["theta"]],
    ds = delta           / spec[["sigma"]],
    db = delta           / spec[["beta"]],
    at = spec[["alpha"]] / spec[["theta"]]
  )
}

# ── Time-domain helpers ───────────────────────────────────────────────────────

.nzc <- function(x) sum(diff(sign(x)) != 0L)

.hjorth <- function(x) {
  d1     <- diff(x)
  d2     <- diff(d1)
  var_x  <- stats::var(x)
  var_d1 <- stats::var(d1)
  var_d2 <- stats::var(d2)
  eps    <- .Machine$double.eps
  hmob   <- sqrt(var_d1 / (var_x  + eps))
  hcomp  <- sqrt(var_d2 / (var_d1 + eps)) / (hmob + eps)
  c(hmob = hmob, hcomp = hcomp)
}

.petrosian_fd <- function(x) {
  N   <- length(x)
  nzc <- .nzc(x)
  log10(N) / (log10(N) + log10(N / (N + 0.4 * nzc)))
}

.higuchi_fd <- function(x, kmax = 10L) {
  N     <- length(x)
  k_seq <- seq_len(kmax)
  Lk    <- vapply(k_seq, function(k) {
    Lmk <- vapply(seq_len(k), function(m) {
      idx <- seq.int(m, N, by = k)
      if (length(idx) < 2L) return(0)
      (N - 1L) / (floor((N - m) / k) * k ^ 2L) * sum(abs(diff(x[idx])))
    }, numeric(1))
    mean(Lmk[Lmk > 0])
  }, numeric(1))
  valid <- Lk > 0
  if (sum(valid) < 2L) return(NA_real_)
  unname(stats::coef(stats::lm(log(Lk[valid]) ~ log(1 / k_seq[valid])))[2L])
}

.perm_entropy <- function(x, order = 3L, delay = 1L) {
  N <- length(x)
  n <- N - (order - 1L) * delay
  if (n <= 0L) return(NA_real_)

  lags     <- seq(0L, (order - 1L) * delay, by = delay)
  idx_mat  <- outer(seq_len(n), lags, "+")
  embedded <- matrix(x[as.vector(idx_mat)], nrow = n, ncol = order)

  pat_str <- apply(t(apply(embedded, 1L, order)), 1L, paste, collapse = "")
  freqs   <- as.vector(table(pat_str)) / n
  h       <- -sum(freqs * log(freqs + .Machine$double.eps))
  h / log(factorial(order))
}

.stat_features <- function(x) {
  mu  <- mean(x)
  sig <- stats::sd(x)
  eps <- .Machine$double.eps
  if (sig < eps) return(c(std = 0, iqr = 0, skew = 0, kurt = 0))
  zx  <- (x - mu) / sig
  c(
    std  = sig,
    iqr  = stats::IQR(x),
    skew = mean(zx ^ 3),
    kurt = mean(zx ^ 4) - 3
  )
}

# ── Per-epoch feature vectors ─────────────────────────────────────────────────
# Accept pre-filtered epoch signals — filtering happens on the full recording
# in .extract_staging_features before epoch extraction (matches YASA structure).

# 21 EEG base features
.eeg_epoch_features <- function(sig, sr) {
  spec   <- .spectral_features(sig, sr)
  ratios <- .spectral_ratios(spec)
  hjorth <- .hjorth(sig)
  stats  <- .stat_features(sig)
  c(
    abspow    = unname(spec["abspow"]),
    alpha     = unname(spec["alpha"]),
    at        = unname(ratios["at"]),
    beta      = unname(spec["beta"]),
    db        = unname(ratios["db"]),
    ds        = unname(ratios["ds"]),
    dt        = unname(ratios["dt"]),
    fdelta    = unname(spec["fdelta"]),
    hcomp     = unname(hjorth["hcomp"]),
    higuchi   = .higuchi_fd(sig),
    hmob      = unname(hjorth["hmob"]),
    iqr       = unname(stats["iqr"]),
    kurt      = unname(stats["kurt"]),
    nzc       = .nzc(sig),
    perm      = .perm_entropy(sig),
    petrosian = .petrosian_fd(sig),
    sdelta    = unname(spec["sdelta"]),
    sigma     = unname(spec["sigma"]),
    skew      = unname(stats["skew"]),
    std       = unname(stats["std"]),
    theta     = unname(spec["theta"])
  )
}

# 17 EOG base features (no ratio features)
.eog_epoch_features <- function(sig, sr) {
  spec   <- .spectral_features(sig, sr)
  hjorth <- .hjorth(sig)
  stats  <- .stat_features(sig)
  c(
    abspow    = unname(spec["abspow"]),
    alpha     = unname(spec["alpha"]),
    beta      = unname(spec["beta"]),
    fdelta    = unname(spec["fdelta"]),
    hcomp     = unname(hjorth["hcomp"]),
    higuchi   = .higuchi_fd(sig),
    hmob      = unname(hjorth["hmob"]),
    iqr       = unname(stats["iqr"]),
    kurt      = unname(stats["kurt"]),
    nzc       = .nzc(sig),
    perm      = .perm_entropy(sig),
    petrosian = .petrosian_fd(sig),
    sdelta    = unname(spec["sdelta"]),
    sigma     = unname(spec["sigma"]),
    skew      = unname(stats["skew"]),
    std       = unname(stats["std"]),
    theta     = unname(spec["theta"])
  )
}

# 11 EMG base features (absolute power + nonlinear only, no spectral bands)
.emg_epoch_features <- function(sig, sr) {
  spec   <- .spectral_features(sig, sr)
  hjorth <- .hjorth(sig)
  stats  <- .stat_features(sig)
  c(
    abspow    = unname(spec["abspow"]),
    hcomp     = unname(hjorth["hcomp"]),
    higuchi   = .higuchi_fd(sig),
    hmob      = unname(hjorth["hmob"]),
    iqr       = unname(stats["iqr"]),
    kurt      = unname(stats["kurt"]),
    nzc       = .nzc(sig),
    perm      = .perm_entropy(sig),
    petrosian = .petrosian_fd(sig),
    skew      = unname(stats["skew"]),
    std       = unname(stats["std"])
  )
}

# ── Normalisation helpers ─────────────────────────────────────────────────────
# YASA normalises via:
#   _c7min_norm : triangular-weighted rolling mean (k=15, centered) → robust_scale
#   _p2min_norm : uniform rolling mean (k=4, right-aligned)         → robust_scale
# robust_scale = (x - median) / (q95 - q5)

.robust_scale <- function(x, q_low = 0.05, q_high = 0.95) {
  eps <- 1e-10
  med <- median(x, na.rm = TRUE)
  q   <- quantile(x, c(q_low, q_high), na.rm = TRUE)
  (x - med) / (q[2L] - q[1L] + eps)
}

# Triangular-weighted rolling mean, centered window of size k.
# Matches pandas rolling(window=k, center=True, min_periods=1, win_type='triang').mean()
.roll_triang_mean <- function(x, k = 15L) {
  n      <- length(x)
  half   <- (k - 1L) %/% 2L
  w_full <- c(seq_len(half + 1L), seq(half, 1L))
  vapply(seq_len(n), function(i) {
    i_start <- max(1L, i - half)
    i_end   <- min(n, i + half)
    w_start <- (i_start - (i - half)) + 1L
    w_end   <- w_start + (i_end - i_start)
    w_sub   <- w_full[w_start:w_end]
    sum(x[i_start:i_end] * w_sub) / sum(w_sub)
  }, numeric(1))
}

# For each base feature column, add _c7min_norm and _p2min_norm, interleaved.
.add_norm_variants <- function(mat, prefix) {
  base_names <- colnames(mat)
  out <- vector("list", ncol(mat) * 3L)
  nms <- character(ncol(mat) * 3L)
  j   <- 1L
  for (i in seq_along(base_names)) {
    nm  <- base_names[i]
    col <- as.vector(mat[, i])

    c7 <- .robust_scale(.roll_triang_mean(col, k = 15L))

    p2_raw <- as.vector(zoo::rollapply(col, 4L, mean, fill = NA,
                                        partial = TRUE, align = "right"))
    p2 <- .robust_scale(p2_raw)

    out[[j]]      <- col; nms[j]      <- paste0(prefix, nm)
    out[[j + 1L]] <- c7;  nms[j + 1L] <- paste0(prefix, nm, "_c7min_norm")
    out[[j + 2L]] <- p2;  nms[j + 2L] <- paste0(prefix, nm, "_p2min_norm")
    j <- j + 3L
  }
  out <- do.call(cbind, out)
  colnames(out) <- nms
  out
}

# NA-fill matrix for a missing channel, preserving the correct column names.
.na_channel_matrix <- function(n, prefix, base_names) {
  nms <- unlist(lapply(base_names, function(nm) {
    c(paste0(prefix, nm),
      paste0(prefix, nm, "_c7min_norm"),
      paste0(prefix, nm, "_p2min_norm"))
  }))
  matrix(NA_real_, nrow = n, ncol = length(nms), dimnames = list(NULL, nms))
}

# ── Main extraction function ──────────────────────────────────────────────────

# Extract the 149-feature staging matrix from an mrpheus_psg object.
# Called internally by stage_epochs(). Not exported.
.extract_staging_features <- function(psg, eeg_ch, eog_ch, emg_ch) {
  cli::cli_alert_info("Extracting staging features ({psg$n_epochs} epochs)...")

  n  <- psg$n_epochs
  sr <- function(ch) psg$channel_map$sample_rate[psg$channel_map$label == ch]

  # Resample to 100 Hz, filter the full recording, then extract epochs.
  # Matches YASA: raw_pick.resample(100) → filter_data() → sliding_window().
  # All features are always computed at 100 Hz regardless of native channel rate.
  SR_TARGET <- 100L

  filter_and_epoch <- function(ch) {
    sig_raw <- psg$edf$signals[[ch]]$signal
    sr_orig <- sr(ch)

    # Resample to 100 Hz if needed — matches YASA's raw_pick.resample(100).
    # Use gsignal::resample (polyphase Kaiser FIR) not stats::approx (linear),
    # to match scipy.signal.resample_poly used internally by MNE.
    if (sr_orig != SR_TARGET) {
      sig_raw <- as.vector(gsignal::resample(sig_raw, SR_TARGET, sr_orig))
    }

    ep_len   <- as.integer(psg$epoch_s * SR_TARGET)
    sig_filt <- .bandpass_filter(sig_raw, SR_TARGET)
    lapply(seq_len(n), function(i) {
      start <- (i - 1L) * ep_len + 1L
      sig_filt[start:(start + ep_len - 1L)]
    })
  }

  # ── EEG ──────────────────────────────────────────────────────────────────
  eeg_epochs <- filter_and_epoch(eeg_ch)
  eeg_mat    <- do.call(rbind, lapply(eeg_epochs, .eeg_epoch_features, sr = SR_TARGET))
  eeg_out    <- .add_norm_variants(eeg_mat, "eeg_")

  # ── Time ─────────────────────────────────────────────────────────────────
  time_hour <- (seq_len(n) - 1L) * psg$epoch_s / 3600
  time_norm <- if (max(time_hour) > 0) time_hour / max(time_hour) else time_hour
  time_out  <- cbind(time_hour = time_hour, time_norm = time_norm)

  # ── EOG ──────────────────────────────────────────────────────────────────
  eog_base <- c("abspow","alpha","beta","fdelta","hcomp","higuchi","hmob",
                "iqr","kurt","nzc","perm","petrosian","sdelta","sigma",
                "skew","std","theta")
  if (!is.na(eog_ch)) {
    eog_epochs <- filter_and_epoch(eog_ch)
    eog_mat    <- do.call(rbind, lapply(eog_epochs, .eog_epoch_features, sr = SR_TARGET))
    eog_out    <- .add_norm_variants(eog_mat, "eog_")
  } else {
    eog_out <- .na_channel_matrix(n, "eog_", eog_base)
  }

  # ── EMG ──────────────────────────────────────────────────────────────────
  emg_base <- c("abspow","hcomp","higuchi","hmob","iqr","kurt","nzc",
                "perm","petrosian","skew","std")
  if (!is.na(emg_ch)) {
    emg_epochs <- filter_and_epoch(emg_ch)
    emg_mat    <- do.call(rbind, lapply(emg_epochs, .emg_epoch_features, sr = SR_TARGET))
    emg_out    <- .add_norm_variants(emg_mat, "emg_")
  } else {
    emg_out <- .na_channel_matrix(n, "emg_", emg_base)
  }

  # ── Assemble and sort alphabetically — matches YASA's features.sort_index() ──
  feat_mat <- cbind(eeg_out, time_out, eog_out, emg_out)
  feat_mat <- feat_mat[, order(colnames(feat_mat))]

  stopifnot(ncol(feat_mat) == 149L)

  tibble::as_tibble(cbind(epoch = seq_len(n), as.data.frame(feat_mat)))
}
