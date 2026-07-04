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
# NOTE: validate against YASA's yasa.SleepStaging._check_params() if staging
# accuracy is lower than expected. Sigma lower bound (12 vs 11 Hz) is the most
# likely source of parity drift.

.STAGING_BANDS <- list(
  sdelta = c(0.5,  2),
  fdelta = c(2,    4),
  theta  = c(4,    8),
  alpha  = c(8,   13),
  sigma  = c(12,  16),
  beta   = c(16,  30)
)

.ABSPOW_RANGE <- c(0.5, 40)

# ── Spectral helpers ──────────────────────────────────────────────────────────

# Compute Welch PSD and return named vector of band features:
# sdelta, fdelta, theta, alpha, sigma, beta (relative) + abspow (absolute).
.spectral_features <- function(sig, sr) {
  nfft   <- 2L ^ ceiling(log2(4 * sr))   # 4-second Welch window
  novlap <- nfft %/% 2L                  # 50 % overlap

  psd <- gsignal::pwelch(sig, fs = sr, window = nfft,
                          noverlap = novlap, nfft = nfft)

  # Absolute total power (0.5–40 Hz) — stored as abspow
  idx_total <- psd$freq >= .ABSPOW_RANGE[1] & psd$freq <= .ABSPOW_RANGE[2]
  abspow    <- pracma::trapz(psd$freq[idx_total], psd$spec[idx_total])

  # Absolute band powers, then normalise to relative
  bp <- vapply(.STAGING_BANDS, function(b) {
    idx <- psd$freq >= b[1] & psd$freq < b[2]
    if (!any(idx)) return(NA_real_)
    pracma::trapz(psd$freq[idx], psd$spec[idx])
  }, numeric(1))

  total_bp <- sum(bp, na.rm = TRUE)
  rel      <- bp / total_bp

  c(rel, abspow = abspow)
}

# Spectral ratio features (EEG only).
.spectral_ratios <- function(spec) {
  c(
    dt = spec[["sdelta"]] / spec[["theta"]],
    ds = spec[["sdelta"]] / spec[["sigma"]],
    db = spec[["sdelta"]] / spec[["beta"]],
    at = spec[["alpha"]]  / spec[["theta"]]
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
  stats::coef(stats::lm(log(Lk[valid]) ~ log(1 / k_seq[valid])))[2L]
}

# Normalised permutation entropy (order = 3, delay = 1, matching YASA defaults)
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
  n   <- length(x)
  mu  <- mean(x)
  sig <- stats::sd(x)
  eps <- .Machine$double.eps
  if (sig < eps) return(c(std = 0, iqr = 0, skew = 0, kurt = 0))
  zx  <- (x - mu) / sig
  c(
    std  = sig,
    iqr  = stats::IQR(x),
    skew = mean(zx ^ 3),
    kurt = mean(zx ^ 4) - 3   # excess kurtosis
  )
}

# ── Per-epoch feature vectors ─────────────────────────────────────────────────
# Base names are alphabetical — must match model feature_names exactly.

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

# ── Rolling normalisation ─────────────────────────────────────────────────────

# Z-score normalise a single feature vector with a rolling window.
# partial = TRUE so edge epochs use whatever data is available (matching pandas
# rolling(min_periods=1) behaviour).
.roll_znorm <- function(x, k, align) {
  eps   <- 1e-10
  rmean <- zoo::rollapply(x, k, mean, fill = NA, partial = TRUE, align = align)
  rsd   <- zoo::rollapply(x, k, function(w) {
    s <- stats::sd(w)
    if (is.na(s) || is.nan(s)) 0 else s
  }, fill = NA, partial = TRUE, align = align)
  rmean[is.na(rmean)] <- x[is.na(rmean)]
  (x - rmean) / (rsd + eps)
}

# For each base feature column, add _c7min_norm (centered 15-epoch window) and
# _p2min_norm (right-aligned 4-epoch window), interleaved: raw, c7min, p2min.
.add_norm_variants <- function(mat, prefix) {
  base_names <- colnames(mat)
  out <- vector("list", ncol(mat) * 3L)
  nms <- character(ncol(mat) * 3L)
  j   <- 1L
  for (i in seq_along(base_names)) {
    nm        <- base_names[i]
    col       <- as.vector(mat[, i])
    out[[j]]     <- col;                         nms[j]     <- paste0(prefix, nm)
    out[[j + 1L]] <- .roll_znorm(col, 15L, "center"); nms[j + 1L] <- paste0(prefix, nm, "_c7min_norm")
    out[[j + 2L]] <- .roll_znorm(col, 4L,  "right");  nms[j + 2L] <- paste0(prefix, nm, "_p2min_norm")
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

  # ── EEG ──────────────────────────────────────────────────────────────────
  eeg_sr  <- sr(eeg_ch)
  eeg_mat <- do.call(rbind, lapply(seq_len(n), function(i)
    .eeg_epoch_features(psg$epochs[[i]][[eeg_ch]], eeg_sr)
  ))
  eeg_out <- .add_norm_variants(eeg_mat, "eeg_")

  # ── Time ─────────────────────────────────────────────────────────────────
  time_hour <- (seq_len(n) - 1L) * psg$epoch_s / 3600
  time_norm <- if (max(time_hour) > 0) time_hour / max(time_hour) else time_hour
  time_out  <- cbind(time_hour = time_hour, time_norm = time_norm)

  # ── EOG ──────────────────────────────────────────────────────────────────
  eog_base <- c("abspow","alpha","beta","fdelta","hcomp","higuchi","hmob",
                "iqr","kurt","nzc","perm","petrosian","sdelta","sigma",
                "skew","std","theta")
  if (!is.na(eog_ch)) {
    eog_sr  <- sr(eog_ch)
    eog_mat <- do.call(rbind, lapply(seq_len(n), function(i)
      .eog_epoch_features(psg$epochs[[i]][[eog_ch]], eog_sr)
    ))
    eog_out <- .add_norm_variants(eog_mat, "eog_")
  } else {
    eog_out <- .na_channel_matrix(n, "eog_", eog_base)
  }

  # ── EMG ──────────────────────────────────────────────────────────────────
  emg_base <- c("abspow","hcomp","higuchi","hmob","iqr","kurt","nzc",
                "perm","petrosian","skew","std")
  if (!is.na(emg_ch)) {
    emg_sr  <- sr(emg_ch)
    emg_mat <- do.call(rbind, lapply(seq_len(n), function(i)
      .emg_epoch_features(psg$epochs[[i]][[emg_ch]], emg_sr)
    ))
    emg_out <- .add_norm_variants(emg_mat, "emg_")
  } else {
    emg_out <- .na_channel_matrix(n, "emg_", emg_base)
  }

  # ── Assemble in model order: EEG (63) | time (2) | EOG (51) | EMG (33) ──
  feat_mat <- cbind(eeg_out, time_out, eog_out, emg_out)

  stopifnot(ncol(feat_mat) == 149L)

  tibble::as_tibble(cbind(epoch = seq_len(n), as.data.frame(feat_mat)))
}
