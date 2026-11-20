#' Detect sleep spindles
#'
#' Detects sleep spindles in an EEG channel with the algorithm of
#' `yasa.spindles_detect()` in YASA (Vallat & Walker, 2021; YASA 0.8.1), which
#' follows Lacourse et al. (2019), and reproduces its output: the same events,
#' at the same sample positions, with the same properties (see Details).
#'
#' The channel is band-passed twice with FIR filters (Hamming window, zero phase,
#' as `mne.filter.filter_data()`) applied to the continuous signal: a broadband
#' filter (`freq_broad`, default 1-30 Hz) and a sigma filter (`freq_spindle`,
#' default 12-15 Hz, transition bands of 1.5 Hz). Three criteria are evaluated at
#' every sample: the relative power of the sigma band in the broadband signal (a
#' short-time Fourier transform with 2 s windows every 0.2 s, linearly
#' interpolated), the moving correlation between the sigma and the broadband
#' signal (0.3 s windows) and the moving RMS of the sigma signal (0.3 s windows;
#' the threshold is the mean plus `thresh_rms` trimmed standard deviations
#' within the accepted epochs, at most 10 uV). The number of criteria met at each
#' sample is smoothed over 0.1 s and a sample belongs to a spindle where at least
#' all but one of the criteria are met. Events closer than `min_distance_ms` are
#' merged, and those whose duration is within `duration_s` are kept.
#'
#' @param psg An `mrpheus_psg` object from [mrpheus::prepare_psg()].
#' @param channel Character. EEG channel label. If `NULL` (default), the first
#'   non-bad EEG channel is used.
#' @param stages Integer vector or `NULL`. Epoch indices in which spindles may
#'   occur (for example the N1, N2 and N3 epochs). If `NULL`, all epochs are
#'   used. As in YASA, the RMS threshold is computed within these epochs and
#'   samples outside them cannot start a spindle, but an event can extend a
#'   little beyond them (the smoothing spreads over 0.1 s).
#' @param freq_spindle Numeric vector of length 2. Spindle band (Hz), as YASA's
#'   `freq_sp`. Default `c(12, 15)`.
#' @param freq_broad Numeric vector of length 2. Broadband filter (Hz). Its
#'   transition bands are chosen as MNE does for `filter_data()`: 25 % of the
#'   cutoff, at least 2 Hz, and no more than the cutoff for the high-pass.
#'   Default `c(1, 30)`.
#' @param duration_s Numeric vector of length 2. Minimum and maximum duration
#'   of a spindle (s), exclusive. Default `c(0.5, 2)`.
#' @param min_distance_ms Numeric or `NULL`. Spindles separated by less than
#'   this (ms) are merged. Default `500`. `NULL` or `0` switches merging off.
#' @param thresh_rel_pow Numeric or `NULL`. Minimum relative sigma power.
#'   Default `0.20`. `NULL` switches the criterion off.
#' @param thresh_corr Numeric or `NULL`. Minimum moving correlation. Default
#'   `0.65`. `NULL` switches the criterion off.
#' @param thresh_rms Numeric or `NULL`. Number of trimmed standard deviations
#'   above the mean for the moving RMS. Default `1.5`. `NULL` switches the
#'   criterion off. At least one criterion must be given.
#'
#' @details
#' Amplitudes must be in uV. As YASA, the channel is skipped (with a warning)
#' when the standard deviation of its central 90 % of samples is outside
#' 0.1-1000 uV, which usually means a wrong unit.
#'
#' The result was compared with `yasa.spindles_detect()` on the Sleep-EDF night
#' SC4001E0 (8.8 h at 100 Hz, and a 6.5 h version at 256 Hz) with a stage mask
#' that cuts through events: the numbers of events are equal (73 and 78), the
#' start and end times are identical, and every property agrees to about 1e-11.
#' The test suite checks the same agreement against stored YASA output on a
#' synthetic EEG with spindles at 100 and 256 Hz
#' (`tests/testthat/test-yasa-spindle-parity.R`).
#'
#' Not implemented: YASA's `multi_only` and `remove_outliers` options, and its
#' spindle - slow oscillation coupling.
#'
#' @return A tibble with one row per detected spindle:
#' \describe{
#'   \item{epoch}{Integer. Epoch of the first sample of the spindle.}
#'   \item{start_s}{Numeric. Onset (s) relative to the start of `epoch`.}
#'   \item{end_s}{Numeric. Offset (s) relative to the start of `epoch`; can
#'     exceed the epoch length when the spindle crosses an epoch boundary.}
#'   \item{duration_s}{Numeric. Duration (s).}
#'   \item{peak_freq_hz}{Numeric. Median instantaneous frequency (Hz) of the
#'     sigma-filtered signal within the spindle (YASA's `Frequency`).}
#'   \item{rms_uv}{Numeric. Root mean square (uV) of the detrended broadband
#'     signal within the spindle (YASA's `RMS`).}
#'   \item{channel}{Character.}
#'   \item{peak_s}{Numeric. Time of the most prominent peak (s from the start of
#'     the first epoch; YASA's `Peak`).}
#'   \item{amplitude_uv}{Numeric. Peak-to-peak amplitude of the detrended
#'     broadband signal (uV; YASA's `Amplitude`).}
#'   \item{amp_filtered_uv}{Numeric. Peak-to-peak amplitude of the sigma-filtered
#'     signal (uV; YASA's `AmpFiltered`).}
#'   \item{abs_power}{Numeric. Median log10 instantaneous power of the sigma
#'     signal (YASA's `AbsPower`).}
#'   \item{rel_power}{Numeric. Median relative sigma power (YASA's `RelPower`).}
#'   \item{oscillations}{Integer. Number of oscillations: peaks at least 60 ms
#'     apart.}
#'   \item{symmetry}{Numeric. Position of the most prominent peak within the
#'     spindle, from 0 to 1 (YASA's `Symmetry`).}
#' }
#'
#' @references
#' Lacourse, K., Delfrate, J., Beaudry, J., Peppard, P., & Warby, S. C. (2019).
#' A sleep spindle detection algorithm that emulates human expert spindle
#' scoring. *Journal of Neuroscience Methods*, 316, 3-11.
#'
#' Vallat, R., & Walker, M. P. (2021). An open-source, high-performance tool
#' for automated sleep staging. *eLife*, 10, e70092.
#' \doi{10.7554/eLife.70092}
#'
#' @export
compute_spindles <- function(psg,
                             channel         = NULL,
                             stages          = NULL,
                             freq_spindle    = c(12, 15),
                             freq_broad      = c(1, 30),
                             duration_s      = c(0.5, 2),
                             min_distance_ms = 500,
                             thresh_rel_pow  = 0.20,
                             thresh_corr     = 0.65,
                             thresh_rms      = 1.5) {
  stopifnot(inherits(psg, "mrpheus_psg"))
  if (is.null(thresh_rel_pow) && is.null(thresh_corr) && is.null(thresh_rms))
    cli::cli_abort("At least one of {.arg thresh_rel_pow}, {.arg thresh_corr} and {.arg thresh_rms} must be given.")

  if (is.null(channel)) {
    channel <- psg$channel_map$label[
      psg$channel_map$type == "EEG" & !psg$channel_map$bad
    ][1]
  }

  sr       <- psg$channel_map$sample_rate[psg$channel_map$label == channel]
  n_epochs <- length(psg$epochs)
  spe      <- as.integer(round(psg$epoch_s * sr))
  epoch_s  <- spe / sr

  # the continuous channel signal, and the sample mask of the accepted epochs
  sig <- unlist(lapply(psg$epochs, function(e) e[[channel]]), use.names = FALSE)
  if (length(sig) != n_epochs * spe)
    cli::cli_abort("Channel {.val {channel}} is not present, with {spe} samples, in every epoch.")
  keep_epoch <- if (is.null(stages)) rep(TRUE, n_epochs) else seq_len(n_epochs) %in% stages
  mask       <- rep(keep_epoch, each = spe)

  empty <- tibble::tibble(
    epoch = integer(), start_s = numeric(), end_s = numeric(), duration_s = numeric(),
    peak_freq_hz = numeric(), rms_uv = numeric(), channel = character(),
    peak_s = numeric(), amplitude_uv = numeric(), amp_filtered_uv = numeric(),
    abs_power = numeric(), rel_power = numeric(), oscillations = integer(),
    symmetry = numeric()
  )

  # YASA skips a channel whose trimmed standard deviation is not 0.1-1000 uV
  trim_sd <- stats::sd(.trim_both(sig, 0.05))
  if (!(trim_sd > 0.1 && trim_sd < 1e3)) {
    cli::cli_warn(c(
      "Wrong data amplitude for channel {.val {channel}} (trimmed SD = {signif(trim_sd, 3)}).",
      "i" = "The unit of the data must be uV. The channel is skipped."
    ))
    return(empty)
  }

  ev <- .yasa_spindles_detect(
    sig, sr, mask, freq_spindle, freq_broad, duration_s, min_distance_ms,
    thresh_rel_pow, thresh_corr, thresh_rms
  )
  if (is.null(ev)) {
    cli::cli_alert_success("Detected 0 spindles in channel {.val {channel}}.")
    return(empty)
  }

  epoch <- as.integer(ev$StartIdx %/% spe) + 1L
  out <- tibble::tibble(
    epoch           = epoch,
    start_s         = ev$Start - (epoch - 1L) * epoch_s,
    end_s           = ev$End   - (epoch - 1L) * epoch_s,
    duration_s      = ev$Duration,
    peak_freq_hz    = ev$Frequency,
    rms_uv          = ev$RMS,
    channel         = channel,
    peak_s          = ev$Peak,
    amplitude_uv    = ev$Amplitude,
    amp_filtered_uv = ev$AmpFiltered,
    abs_power       = ev$AbsPower,
    rel_power       = ev$RelPower,
    oscillations    = as.integer(ev$Oscillations),
    symmetry        = ev$Symmetry
  )
  cli::cli_alert_success(
    "Detected {nrow(out)} spindles in channel {.val {channel}}."
  )
  out
}

# ── Internals: port of yasa.spindles_detect() (YASA 0.8.1) ────────────────────

# MNE band-pass with automatic transition bands:
# l_trans = min(max(0.25 * l_freq, 2), l_freq), h_trans = min(max(0.25 * h_freq, 2), sfreq / 2 - h_freq)
.mne_auto_trans <- function(sfreq, l_freq, h_freq) {
  c(l = min(max(0.25 * l_freq, 2), l_freq), h = min(max(0.25 * h_freq, 2), sfreq / 2 - h_freq))
}

# scipy.fftpack.next_fast_len: the next number whose only prime factors are 2, 3 and 5
.next_fast_len <- function(n) {
  repeat {
    m <- n
    for (p in c(2, 3, 5)) while (m %% p == 0) m <- m %/% p
    if (m == 1) return(as.integer(n))
    n <- n + 1
  }
}

# scipy.signal.hilbert(x, N)[:length(x)]
.hilbert_analytic <- function(x, N) {
  n0 <- length(x)
  X  <- stats::fft(c(x, numeric(N - n0)))
  h  <- numeric(N)
  if (N %% 2L == 0L) { h[1L] <- 1; h[N / 2L + 1L] <- 1; h[2L:(N / 2L)] <- 2 }
  else               { h[1L] <- 1; h[2L:((N + 1L) / 2L)] <- 2 }
  stats::fft(X * h, inverse = TRUE)[seq_len(n0)] / N
}

# np.interp / scipy interp1d(kind = "linear", fill_value = 0) at `x`
.lin_interp0 <- function(xp, fp, x) {
  j   <- findInterval(x, xp)                    # xp[j] <= x < xp[j + 1]
  out <- numeric(length(x))
  n   <- length(xp)
  inside <- x >= xp[1L] & x <= xp[n]
  j[j == n] <- n - 1L                           # x == xp[n]
  k <- which(inside & j >= 1L)
  slope <- (fp[j[k] + 1L] - fp[j[k]]) / (xp[j[k] + 1L] - xp[j[k]])
  out[k] <- slope * (x[k] - xp[j[k]]) + fp[j[k]]
  out
}

# yasa.stft_power(x, sf, window = 2, step = 0.2, band, interp = FALSE): relative
# power of the sigma band (freq_sp) among `band`, one value per STFT frame
.stft_rel_pow <- function(x, sf, freq_sp, band, window = 2, step = 0.2) {
  nperseg  <- as.integer(window * sf)
  noverlap <- as.integer(nperseg - (step * sf))
  nstep    <- nperseg - noverlap
  # scipy.signal.stft: zero boundary of nperseg %/% 2, then zero padding to a whole number of steps
  half <- nperseg %/% 2L
  N0   <- length(x) + 2L * half
  nadd <- ((-(N0 - nperseg)) %% nstep) %% nperseg
  L    <- N0 + nadd
  nseg <- (L - noverlap) %/% nstep
  # frame times and frequencies, with scipy's own float arithmetic
  tt <- (nperseg / 2 + (0:(ceiling((L - nperseg / 2 + 1 - nperseg / 2) / nstep) - 1)) * nstep) / sf - (nperseg / 2) / sf
  stopifnot(length(tt) == nseg)
  d  <- 1 / sf
  f  <- (0:(nperseg %/% 2L)) * (1 / (nperseg * d))
  in_band  <- f >= band[1L] & f <= band[2L]
  in_sigma <- f[in_band] >= freq_sp[1L] & f[in_band] <= freq_sp[2L]
  win <- 0.5 - 0.5 * cos(2 * pi * (0:(nperseg - 1L)) / nperseg)           # periodic Hann
  xe  <- c(numeric(half), x, numeric(L - length(x) - half))
  rel <- numeric(nseg)
  chunk <- 20000L
  for (s in seq.int(1L, nseg, by = chunk)) {
    e   <- min(s + chunk - 1L, nseg)
    idx <- outer(seq_len(nperseg), (s - 1L):(e - 1L) * nstep, `+`)
    M   <- matrix(xe[idx], nrow = nperseg) * win
    P   <- Mod(stats::mvfft(M)[which(in_band), , drop = FALSE])^2
    rel[s:e] <- colSums(P[in_sigma, , drop = FALSE]) / colSums(P)
  }
  list(t = tt, rel = rel)
}

# yasa.moving_transform(x, y, sf, window = 0.3, step = 0.1, method, interp = TRUE)
.moving_transform <- function(x, y = NULL, sf, method = c("corr", "rms"), window = 0.3, step = 0.1) {
  method  <- match.arg(method)
  n       <- length(x)
  halfdur <- window / 2
  total   <- n / sf
  last    <- n - 1L
  idx     <- (0:(ceiling((total - 0) / step) - 1)) * step
  beg     <- as.integer(trunc((idx - halfdur) * sf))
  end     <- as.integer(trunc((idx + halfdur) * sf))
  beg[beg < 0L] <- 0L
  end[end > last] <- n
  t   <- (beg + end) / 2 / sf
  win <- as.numeric(end - beg)
  cs  <- function(v) c(0, cumsum(v))
  if (method == "rms") {
    cx2 <- cs(x^2)
    out <- ifelse(win > 0, sqrt((cx2[end + 1L] - cx2[beg + 1L]) / win), NA_real_)
  } else {
    cx <- cs(x); cy <- cs(y); cxy <- cs(x * y); cx2 <- cs(x^2); cy2 <- cs(y^2)
    sx  <- cx[end + 1L]  - cx[beg + 1L];  sy  <- cy[end + 1L]  - cy[beg + 1L]
    sxy <- cxy[end + 1L] - cxy[beg + 1L]
    sx2 <- cx2[end + 1L] - cx2[beg + 1L]; sy2 <- cy2[end + 1L] - cy2[beg + 1L]
    num <- sxy - sx * sy / win
    vx  <- pmax(sx2 - sx^2 / win, 0); vy <- pmax(sy2 - sy^2 / win, 0)
    den <- sqrt(vx * vy)
    out <- ifelse(win >= 2 & den != 0, num / den, NA_real_)
  }
  .lin_interp0(t, out, (0:(n - 1L)) / sf)
}

# The moving-RMS threshold of yasa.spindles_detect(): mean plus `mult` trimmed
# standard deviations (10 % cut at each end) of the values in the stage mask,
# at most 10 uV (to avoid too high a threshold caused by artefacts or movement)
.sp_rms_threshold <- function(mrms, mask, mult) {
  v  <- sort(mrms[mask])
  lc <- as.integer(0.10 * length(v))
  min(mean(mrms[mask]) + mult * stats::sd(v[seq.int(lc + 1L, length(v) - lc)]), 10)
}

# yasa.detection._merge_close
.merge_close <- function(index, min_distance_ms, sf) {
  min_distance <- min_distance_ms / 1000 * sf
  d   <- diff(index)
  gap <- which(d > 1L & d < min_distance)
  if (!length(gap)) return(index)
  fill <- unlist(lapply(gap, function(j) seq.int(index[j] + 1L, index[j + 1L] - 1L)))
  sort(c(index, fill))
}

# yasa.detection._detrend_linear
.detrend_linear <- function(x) {
  t <- seq_along(x) - 1 - (length(x) - 1) / 2
  x - mean(x) - t * sum(t * x) / sum(t * t)
}

# scipy.signal.find_peaks(x, distance = distance, prominence = (None, None)):
# 1-based indices of the peaks that survive the distance rule, and their prominences
.peaks_prominence <- function(x, distance) {
  pk <- .sw_find_peaks(x) + 1L                          # 1-based local maxima
  if (!length(pk)) return(list(peaks = integer(0), prom = numeric(0)))
  if (length(pk) > 1L) {                                # scipy _select_by_peak_distance
    dist  <- ceiling(distance)
    keep  <- rep(TRUE, length(pk))
    order <- order(x[pk])                               # ascending priority
    for (i in rev(seq_along(order))) {
      j <- order[i]
      if (!keep[j]) next
      k <- j - 1L
      while (k >= 1L && pk[j] - pk[k] < dist) { keep[k] <- FALSE; k <- k - 1L }
      k <- j + 1L
      while (k <= length(pk) && pk[k] - pk[j] < dist) { keep[k] <- FALSE; k <- k + 1L }
    }
    pk <- pk[keep]
  }
  n <- length(x)
  prom <- vapply(pk, function(p) {
    lmin <- x[p]; i <- p
    while (i >= 1L && x[i] <= x[p]) { if (x[i] < lmin) lmin <- x[i]; i <- i - 1L }
    rmin <- x[p]; i <- p
    while (i <= n && x[i] <= x[p]) { if (x[i] < rmin) rmin <- x[i]; i <- i + 1L }
    x[p] - max(lmin, rmin)
  }, numeric(1))
  list(peaks = pk, prom = prom)
}

# One spindle criterion as 0/1 per sample. A non-finite value (NaN or NA, e.g. a
# moving correlation over a flat stretch) counts as 0, as in YASA's
# `(x >= thresh).astype(int)`; left as NA it would turn the cumulative sum below
# into NA for every later sample.
.sp_crit <- function(v, thr) {
  r <- v >= thr
  r[is.na(r)] <- FALSE
  as.numeric(r)
}

.yasa_spindles_detect <- function(x, sf, mask, freq_sp, freq_broad, duration, min_distance,
                                  thresh_rel_pow, thresh_corr, thresh_rms) {
  n <- length(x)
  tb <- .mne_auto_trans(sf, freq_broad[1L], freq_broad[2L])
  xb <- .mne_fir_apply(x, .mne_fir_bandpass_design(sf, freq_broad[1L], freq_broad[2L], tb[["l"]], tb[["h"]]))
  xs <- .mne_fir_apply(x, .mne_fir_bandpass_design(sf, freq_sp[1L], freq_sp[2L], 1.5, 1.5))

  nfast    <- .next_fast_len(n)
  analytic <- .hilbert_analytic(xs, nfast)
  inst_pow  <- Re(analytic)^2 + Im(analytic)^2
  inst_freq <- sf / (2 * pi) * diff(Arg(analytic))
  rm(analytic)

  st <- .stft_rel_pow(xb, sf, freq_sp, freq_broad)
  do_rel <- !is.null(thresh_rel_pow); do_corr <- !is.null(thresh_corr); do_rms <- !is.null(thresh_rms)
  n_thresh <- do_rel + do_corr + do_rms
  idx_sum <- numeric(n)
  if (do_rel)  idx_sum <- idx_sum + .sp_crit(.lin_interp0(st$t, st$rel, (0:(n - 1L)) / sf), thresh_rel_pow)
  if (do_corr) idx_sum <- idx_sum + .sp_crit(.moving_transform(xs, xb, sf, "corr"), thresh_corr)
  if (do_rms) {
    mrms <- .moving_transform(xs, NULL, sf, "rms")
    idx_sum <- idx_sum + .sp_crit(mrms, .sp_rms_threshold(mrms, mask, thresh_rms))
  }
  idx_sum[!mask] <- 0

  w  <- as.integer(0.1 * sf)                         # np.convolve(idx_sum, ones(w), "same") / w
  cs <- c(0, cumsum(idx_sum))
  lo <- pmax(seq_len(n) - (w %/% 2L), 1L)
  hi <- pmin(seq_len(n) + ((w - 1L) %/% 2L), n)
  sm <- (cs[hi + 1L] - cs[lo]) / w
  where_sp <- which(sm > (n_thresh - 1)) - 1L        # 0-based sample indices
  if (!length(where_sp)) return(NULL)
  if (!is.null(min_distance) && min_distance > 0) where_sp <- .merge_close(where_sp, min_distance, sf)

  grp   <- cumsum(c(1L, diff(where_sp) != 1L))
  first <- where_sp[!duplicated(grp)]
  lastv <- where_sp[!duplicated(grp, fromLast = TRUE)]
  sp_start <- first / sf; sp_end <- lastv / sf; sp_dur <- sp_end - sp_start
  good <- sp_dur > duration[1L] & sp_dur < duration[2L]
  if (!any(good)) return(NULL)
  first <- first[good]; lastv <- lastv[good]
  sp_start <- sp_start[good]; sp_end <- sp_end[good]; sp_dur <- sp_dur[good]

  distance <- 60 * sf / 1000
  m <- length(first)
  res <- data.frame(Start = sp_start, Peak = NA_real_, End = sp_end, Duration = sp_dur,
                    Amplitude = NA_real_, AmpFiltered = NA_real_, RMS = NA_real_,
                    AbsPower = NA_real_, RelPower = NA_real_, Frequency = NA_real_,
                    Oscillations = NA_real_, Symmetry = NA_real_, StartIdx = first)
  for (j in seq_len(m)) {
    ii  <- (first[j] + 1L):(lastv[j] + 1L)           # 1-based
    det <- .detrend_linear(xb[ii])
    res$Amplitude[j]   <- diff(range(det))
    res$AmpFiltered[j] <- diff(range(xs[ii]))
    res$RMS[j]         <- sqrt(mean(det^2))
    i0 <- sum(st$t < sp_start[j]) + 1L               # searchsorted left (1-based start)
    i1 <- sum(st$t <= sp_end[j])                     # searchsorted right, exclusive end -> last 1-based
    res$RelPower[j] <- if (i0 <= i1) stats::median(st$rel[i0:i1]) else st$rel[which.min(abs(st$t - 0.5 * (sp_start[j] + sp_end[j])))]
    ip <- inst_pow[ii]; if_ <- inst_freq[ii]
    res$AbsPower[j]  <- stats::median(log10(ip[ip > 0]))
    res$Frequency[j] <- stats::median(if_[if_ > 0 & !is.na(if_)])
    pp <- .peaks_prominence(det, distance)
    res$Oscillations[j] <- length(pp$peaks)
    if (length(pp$peaks)) {
      pk <- pp$peaks[which.max(pp$prom)] - 1L        # 0-based sample within the spindle
      res$Peak[j]     <- sp_start[j] + pk / sf
      res$Symmetry[j] <- pk / length(det)
    }
  }
  res
}
