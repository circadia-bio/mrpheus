# R/staging_features.R
#
# Internal feature extraction for AASM sleep staging.
# All functions are package-private (prefixed with `.`) — not exported.
#
# Reproduces the YASA feature pipeline in pure R using base R, gsignal,
# and pracma. The 149-feature matrix must match the Python implementation
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
  if (n == 0L) return(x)
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

# ── MNE-equivalent FFT resampling ────────────────────────────────────────────────────────
# Matches mne.filter.resample(x, up=up, down=1, npad='auto', method='fft', pad='auto').
# Called by YASA as raw_pick.resample(100, npad='auto').
#
# Algorithm (from MNE _resample_fft + _smart_pad):
#  1. npad='auto': pad signal to the next power of 2 (makes both FFTs fast)
#  2. reflect_limited padding at each end
#  3. FFT; scale all bins by new_len/orig_len (boxcar window = no windowing)
#  4. Halve the Nyquist bin (upsampling, shorter=FALSE)
#  5. Zero-pad in the frequency domain to new_len (= ideal sinc interpolation)
#  6. IFFT then trim the padding back off

# MNE _smart_pad with pad='reflect_limited':
#   Odd reflection limited to the signal length, zero-padded if npad > len(x)-1.
.mne_smart_pad <- function(x, n_pad) {
  nx    <- length(x)
  l     <- n_pad[1L]
  r     <- n_pad[2L]
  l_lim <- min(l, nx - 1L)
  r_lim <- min(r, nx - 1L)
  c(
    rep(0.0, max(l - nx + 1L, 0L)),             # left zero-pad (rarely needed)
    if (l_lim > 0L) 2.0 * x[1L] - x[seq.int(l_lim + 1L, 2L, by = -1L)] else NULL,
    x,
    if (r_lim > 0L) 2.0 * x[nx] - x[seq.int(nx - 1L, nx - r_lim, by = -1L)] else NULL,
    rep(0.0, max(r - nx + 1L, 0L))              # right zero-pad
  )
}

.mne_fft_resample <- function(x, up) {
  ratio     <- as.double(up)
  x_len     <- length(x)
  final_len <- as.integer(round(x_len * ratio))

  # npad='auto': pad to next power of 2 for fast FFTs
  min_add    <- min(x_len %/% 8L, 100L) * 2L
  n_total    <- 2L ^ ceiling(log2(x_len + min_add))
  npad_total <- n_total - x_len
  npad_left  <- npad_total %/% 2L
  npad_right <- npad_total - npad_left
  npads      <- c(npad_left, npad_right)

  # reflect_limited padding -> power-of-2 length (fast FFT)
  x_pad    <- .mne_smart_pad(x, npads)
  orig_len <- length(x_pad)                      # = n_total
  new_len  <- as.integer(max(round(ratio * orig_len), 1L))

  to_rm_start <- as.integer(round(ratio * npads[1L]))
  to_rm_end   <- new_len - final_len - to_rm_start

  # Forward FFT of padded signal (power-of-2 -> fast)
  X <- stats::fft(x_pad)

  # Boxcar window: scale uniformly by new_len/orig_len
  X <- X * (new_len / orig_len)

  # Nyquist bin halved for upsampling (shorter=FALSE, use_len=orig_len, even length)
  if (orig_len %% 2L == 0L)
    X[orig_len %/% 2L + 1L] <- X[orig_len %/% 2L + 1L] * 0.5

  # Zero-pad in frequency domain (ideal sinc upsampling)
  n_pos        <- orig_len %/% 2L + 1L
  n_neg_total  <- orig_len - n_pos + 1L
  X_new <- complex(length.out = new_len)
  X_new[seq_len(n_pos)] <- X[seq_len(n_pos)]
  if (n_neg_total > 0L)
    X_new[seq.int(new_len - n_neg_total + 1L, new_len)] <- X[seq.int(n_pos, orig_len)]

  # Inverse FFT then trim
  y_full <- Re(stats::fft(X_new, inverse = TRUE)) / new_len
  y_full[seq.int(to_rm_start + 1L, new_len - to_rm_end)]
}

# ── Spectral helpers ──────────────────────────────────────────────────────────

# Composite Simpson's rule matching scipy.integrate.simpson.
# - Odd  N: standard 1/3 rule.
# - Even N: scipy 1.11+ approach — standard 1/3 on first N-3 points,
#            Simpson's 3/8 rule on the last 4 points.
.scipy_simpson <- function(y, dx) {
  N <- length(y)
  if (N < 2L) return(0)
  if (N == 2L) return(dx * (y[1L] + y[2L]) / 2)

  if (N %% 2L == 1L) {
    coef      <- rep(2, N)
    coef[seq(2L, N - 1L, by = 2L)] <- 4
    coef[c(1L, N)] <- 1
    return(dx / 3 * sum(coef * y))
  } else {
    n_std   <- N - 3L
    s_total <- 0
    if (n_std >= 3L) {
      coef <- rep(2, n_std)
      coef[seq(2L, n_std - 1L, by = 2L)] <- 4
      coef[c(1L, n_std)] <- 1
      s_total <- dx / 3 * sum(coef * y[seq_len(n_std)])
    }
    s_total + 3 * dx / 8 * (y[N - 3L] + 3 * y[N - 2L] + 3 * y[N - 1L] + y[N])
  }
}

# Bias correction factor matching scipy.signal.spectral._median_bias(n).
.median_bias <- function(n) {
  n_half <- (n - 1L) %/% 2L
  if (n_half < 1L) return(1)
  ii_2 <- 2 * seq_len(n_half)
  1 + sum(1 / (ii_2 + 1) - 1 / ii_2)
}

# Compute Welch PSD with bias-corrected median averaging over segments.
# Matches YASA/scipy.signal.welch(average='median').
# Uses mvfft() for a single batched C-level FFT over all segments;
# rowmedian_cpp replaces apply(pgrams, 1L, stats::median).
.welch_median_psd <- function(sig, sr, win_n, overlap = 0.5) {
  # Periodic (DFT-even) Hamming window — matches scipy.signal.get_window('hamming', N, fftbins=True)
  wvec   <- 0.54 - 0.46 * cos(2 * pi * seq.int(0L, win_n - 1L) / win_n)
  hop    <- max(1L, as.integer(win_n * (1 - overlap)))
  n_freq <- win_n %/% 2L + 1L
  freqs  <- seq(0, sr / 2, length.out = n_freq)
  if (length(sig) < win_n)
    return(list(freq = freqs, spec = rep(NA_real_, n_freq)))
  starts <- seq(1L, length(sig) - win_n + 1L, by = hop)
  if (length(starts) == 0L)
    return(list(freq = freqs, spec = rep(NA_real_, n_freq)))

  scale  <- sum(wvec ^ 2) * sr
  n_segs <- length(starts)

  idx_mat <- matrix(
    rep(seq_len(win_n), n_segs) + rep((starts - 1L), each = win_n),
    nrow = win_n
  )
  seg_mat <- matrix(sig[idx_mat], nrow = win_n)
  seg_mat <- sweep(seg_mat, 2L, colMeans(seg_mat), "-") * wvec

  ft_mat  <- mvfft(seg_mat)[seq_len(n_freq), , drop = FALSE]
  pgrams  <- Mod(ft_mat) ^ 2 / scale
  pgrams[seq(2L, n_freq - 1L), ] <- 2 * pgrams[seq(2L, n_freq - 1L), ]

  bias    <- .median_bias(n_segs)
  psd_med <- rowmedian_cpp(pgrams) / bias
  list(freq = freqs, spec = psd_med)
}

# Compute Welch PSD and return named vector of band features.
.spectral_features <- function(sig, sr) {
  win <- as.integer(5L * sr)
  psd <- .welch_median_psd(sig, sr, win_n = win)
  dx  <- psd$freq[2L] - psd$freq[1L]

  bp <- vapply(.STAGING_BANDS, function(b) {
    idx <- psd$freq >= b[1] & psd$freq <= b[2]
    if (!any(idx)) return(NA_real_)
    .scipy_simpson(psd$spec[idx], dx)
  }, numeric(1))

  idx_broad_full <- psd$freq >= .FREQ_BROAD[1] & psd$freq <= .FREQ_BROAD[2]
  total_power    <- .scipy_simpson(psd$spec[idx_broad_full], dx)
  rel            <- bp / total_power

  abspow <- pracma::trapz(psd$freq[idx_broad_full], psd$spec[idx_broad_full])
  c(rel, abspow = abspow)
}

# Spectral ratio features (EEG only). Numerator uses full delta (0.4-4 Hz).
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

.nzc <- function(x) nzc_cpp(x)

.hjorth <- function(x) hjorth_cpp(x)

.petrosian_fd <- function(x) petrosian_fd_cpp(x)

.higuchi_fd <- function(x, kmax = 10L) higuchi_fd_cpp(x, kmax)

.perm_entropy <- function(x, order = 3L, delay = 1L) perm_entropy_cpp(x, order, delay)

.roll_triang_mean <- function(x, k = 15L) roll_triang_mean_cpp(x, k)

.stat_features <- function(x) stat_features_cpp(x)

# ── Per-epoch feature vectors ─────────────────────────────────────────────────

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

# 11 EMG base features (absolute power + nonlinear only)
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
#   _c7min_norm : triangular-weighted rolling mean (k=15, centered) -> robust_scale
#   _p2min_norm : uniform rolling mean (k=4, right-aligned)         -> robust_scale
# robust_scale = (x - median) / (q95 - q5)

.robust_scale <- function(x, q_low = 0.05, q_high = 0.95) {
  robust_scale_cpp(x, q_low, q_high)
}

# .roll_triang_mean is defined in the time-domain helpers section above
# (calls roll_triang_mean_cpp). Pure-R reference removed to avoid shadowing.

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

    p2_raw <- roll_right_mean_cpp(col, 4L)
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

  SR_TARGET <- 100L

  filter_and_epoch <- function(ch) {
    sig_raw <- psg$edf$signals[[ch]]$signal
    sr_orig <- sr(ch)
    if (sr_orig != SR_TARGET) {
      sig_raw <- .mne_fft_resample(sig_raw, up = SR_TARGET / sr_orig)
    }
    ep_len   <- as.integer(psg$epoch_s * SR_TARGET)
    sig_filt <- .bandpass_filter(sig_raw, SR_TARGET)
    lapply(seq_len(n), function(i) {
      start <- (i - 1L) * ep_len + 1L
      sig_filt[start:(start + ep_len - 1L)]
    })
  }

  # EEG
  eeg_epochs <- filter_and_epoch(eeg_ch)
  eeg_mat    <- do.call(rbind, lapply(eeg_epochs, .eeg_epoch_features, sr = SR_TARGET))
  eeg_out    <- .add_norm_variants(eeg_mat, "eeg_")

  # Time
  time_hour <- (seq_len(n) - 1L) * psg$epoch_s / 3600
  time_norm <- if (max(time_hour) > 0) time_hour / max(time_hour) else time_hour
  time_out  <- cbind(time_hour = time_hour, time_norm = time_norm)

  # EOG
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

  # EMG
  emg_base <- c("abspow","hcomp","higuchi","hmob","iqr","kurt","nzc",
                "perm","petrosian","skew","std")
  if (!is.na(emg_ch)) {
    emg_epochs <- filter_and_epoch(emg_ch)
    emg_mat    <- do.call(rbind, lapply(emg_epochs, .emg_epoch_features, sr = SR_TARGET))
    emg_out    <- .add_norm_variants(emg_mat, "emg_")
  } else {
    emg_out <- .na_channel_matrix(n, "emg_", emg_base)
  }

  # Assemble and sort alphabetically — matches YASA's features.sort_index()
  feat_mat <- cbind(eeg_out, time_out, eog_out, emg_out)
  feat_mat <- feat_mat[, order(colnames(feat_mat))]

  stopifnot(ncol(feat_mat) == 149L)

  tibble::as_tibble(cbind(epoch = seq_len(n), as.data.frame(feat_mat)))
}
