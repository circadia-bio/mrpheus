# R/compute_band_power.R
#
# Band power per epoch and channel using a Welch PSD that matches
# scipy.signal.welch with YASA's default parameters:
#   win_sec=4, noverlap=200, nfft=1024, window='hann',
#   detrend='constant', average='mean'.
#
# Reference: yasa.bandpower() — Vallat & Walker (2021), eLife 10:e70092.

# ── Internal Welch PSD ────────────────────────────────────────────────────────
#
# Matches scipy.signal.welch(x, fs, window='hann', nperseg, noverlap, nfft,
#   detrend='constant', average='mean', scaling='density').
# Returns list(freq, spec) — one-sided PSD in V^2/Hz.
.welch_psd_bp <- function(x, fs, win_sec, noverlap, nfft) {
  nperseg <- as.integer(win_sec * fs)
  nperseg <- min(nperseg, length(x))

  # Raise nfft to next power of 2 above nperseg if needed (scipy raises error;
  # we silently correct to avoid crashing on unusual sample-rate combinations).
  if (nfft < nperseg) nfft <- 2L ^ ceiling(log2(nperseg))

  # Hann window: symmetric, matches scipy.signal.windows.hann(nperseg).
  k   <- seq(0L, nperseg - 1L)
  win <- 0.5 * (1 - cos(2 * pi * k / (nperseg - 1L)))

  step   <- max(nperseg - as.integer(noverlap), 1L)
  n      <- length(x)
  starts <- seq(1L, max(1L, n - nperseg + 1L), by = step)

  # scipy density scaling: scale = fs * sum(w^2)
  scale <- fs * sum(win^2)
  n_one <- nfft %/% 2L + 1L

  # Compute one-sided periodogram for each segment
  pgrams <- vapply(starts, function(s) {
    seg   <- x[s:(s + nperseg - 1L)]
    seg   <- seg - mean(seg)                         # detrend='constant'
    seg_w <- c(seg * win, numeric(nfft - nperseg))   # apply window, zero-pad
    X     <- stats::fft(seg_w)
    P     <- (Mod(X)^2) / scale

    # One-sided: double all interior bins.
    # Even nfft: DC (index 1) and Nyquist (index n_one) are NOT doubled.
    # Odd nfft:  DC (index 1) is not doubled; no exact Nyquist bin.
    P_one <- P[seq_len(n_one)]
    if (nfft %% 2L == 0L) {
      if (n_one > 2L) P_one[2L:(n_one - 1L)] <- 2 * P_one[2L:(n_one - 1L)]
    } else {
      if (n_one > 1L) P_one[2L:n_one] <- 2 * P_one[2L:n_one]
    }
    P_one
  }, numeric(n_one))

  # Mean across segments (vapply returns n_one x n_segs matrix)
  psd <- if (is.matrix(pgrams)) rowMeans(pgrams) else pgrams

  list(
    freq = seq(0L, n_one - 1L) * (fs / nfft),
    spec = psd
  )
}

# ── Public function ───────────────────────────────────────────────────────────

#' Compute EEG band power per epoch
#'
#' Estimates power spectral density (PSD) using Welch's method and integrates
#' power within standard EEG frequency bands for each epoch and channel.
#'
#' The Welch implementation matches `scipy.signal.welch` with YASA's default
#' parameters (`win_sec = 4`, `noverlap = 200`, `nfft = 1024`, Hann window,
#' constant detrend, mean averaging). Relative power matches YASA's
#' `relative = True` behaviour: each band is divided by the sum of all band
#' powers (not the total PSD integral).
#'
#' @param psg An `mrpheus_psg` object from [mrpheus::prepare_psg()].
#' @param channels Character vector. Channel labels to include. If `NULL`
#'   (default), all non-bad EEG channels are used.
#' @param bands Named list of length-2 numeric vectors defining frequency bands
#'   in Hz. Default matches YASA's `bandpower()`:
#'   ```
#'   list(delta = c(0.5, 4), theta = c(4, 8), alpha = c(8, 12),
#'        sigma = c(12, 16), beta = c(16, 30), gamma = c(30, 45))
#'   ```
#' @param relative Logical. If `TRUE`, each band power is divided by the sum
#'   of all band powers (dimensionless). Matches YASA `relative = True`.
#'   Default `FALSE`.
#' @param win_sec Numeric. Welch window length in seconds. Default `4`.
#' @param noverlap Integer. Number of overlapping samples between Welch windows.
#'   Default `200`. Must be less than `win_sec * sample_rate`.
#' @param nfft Integer. FFT length. Default `1024`. If smaller than the window
#'   length in samples it is automatically raised to the next power of 2.
#'
#' @return A tibble with columns `epoch`, `channel`, one column per named band,
#'   and `total_power` (sum of all band powers, in V^2/Hz units before any
#'   relative scaling).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' bp <- compute_band_power(psg)
#' bp <- compute_band_power(psg, channels = "EEG Fpz-Cz", relative = TRUE)
#'
#' # Custom bands
#' bp <- compute_band_power(
#'   psg,
#'   bands = list(slow = c(0.5, 1), delta = c(1, 4), theta = c(4, 8))
#' )
#' }
compute_band_power <- function(psg,
                               channels = NULL,
                               bands    = list(
                                 delta = c(0.5,  4),
                                 theta = c(4,    8),
                                 alpha = c(8,   12),
                                 sigma = c(12,  16),
                                 beta  = c(16,  30),
                                 gamma = c(30,  45)
                               ),
                               relative = FALSE,
                               win_sec  = 4,
                               noverlap = 200L,
                               nfft     = 1024L) {
  stopifnot(inherits(psg, "mrpheus_psg"))

  if (is.null(channels)) {
    channels <- psg$channel_map$label[
      psg$channel_map$type == "EEG" & !psg$channel_map$bad
    ]
  }

  rows <- lapply(seq_along(psg$epochs), function(i) {
    ep <- psg$epochs[[i]]

    lapply(channels, function(ch) {
      sig <- ep[[ch]]
      sr  <- psg$channel_map$sample_rate[psg$channel_map$label == ch]
      if (is.null(sig) || length(sig) < 2L) return(NULL)

      psd <- .welch_psd_bp(sig, fs = sr,
                            win_sec  = win_sec,
                            noverlap = noverlap,
                            nfft     = nfft)

      band_power <- vapply(bands, function(b) {
        idx <- psd$freq >= b[1] & psd$freq <= b[2]
        if (!any(idx)) return(NA_real_)
        pracma::trapz(psd$freq[idx], psd$spec[idx])
      }, numeric(1))

      total <- sum(band_power, na.rm = TRUE)
      if (relative) band_power <- band_power / total

      row              <- as.list(band_power)
      row$epoch        <- i
      row$channel      <- ch
      row$total_power  <- total
      tibble::as_tibble(row)
    })
  })

  dplyr::bind_rows(unlist(rows, recursive = FALSE)) |>
    dplyr::relocate(epoch, channel)
}
