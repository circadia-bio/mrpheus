# R/compute_band_power.R
#
# Band power per epoch and channel, following yasa.bandpower() (YASA 0.8.1):
# scipy.signal.welch(nperseg = int(win_sec * fs), average = "median",
# window = "hamming") with scipy's own defaults for the rest (overlap of half a
# window, nfft = nperseg, constant detrend, density scaling), Simpson
# integration, and relative power taken as the share of the PSD integral
# between the lowest and highest band edge.
#
# Reference: yasa.bandpower() - Vallat & Walker (2021), eLife 10:e70092.

# ── Internal Welch PSD ────────────────────────────────────────────────────────
#
# Matches scipy.signal.welch(x, fs, window, nperseg, noverlap, nfft,
#   detrend = 'constant', average, scaling = 'density').
#
# `window`   "hann" or "hamming".
# `periodic` TRUE builds the window as scipy.signal.welch does (the periodic,
#            "DFT-even" form, denominator N). FALSE gives the symmetric form
#            (denominator N - 1), which is what the helper produced before
#            compute_band_power() was aligned with YASA; compute_temporal_bandpower()
#            still uses it.
# `average`  "mean", or "median" with scipy's bias correction (.median_bias()).
# `noverlap` NULL = half a window (scipy's default); `nfft` NULL = window length.
#
# Returns list(freq, spec) - one-sided PSD in V^2/Hz.
.welch_psd_bp <- function(x, fs, win_sec, noverlap = NULL, nfft = NULL,
                          window = "hann", average = "mean", periodic = FALSE) {
  window  <- match.arg(window,  c("hann", "hamming"))
  average <- match.arg(average, c("mean", "median"))

  nperseg <- as.integer(win_sec * fs)
  nperseg <- min(nperseg, length(x))

  if (is.null(noverlap)) noverlap <- nperseg %/% 2L
  if (is.null(nfft))     nfft     <- nperseg

  # Raise nfft to next power of 2 above nperseg if needed (scipy raises error;
  # we silently correct to avoid crashing on unusual sample-rate combinations).
  if (nfft < nperseg) nfft <- 2L ^ ceiling(log2(nperseg))

  # Window. Periodic form matches scipy.signal.get_window(name, nperseg) as used
  # inside scipy.signal.welch; symmetric form matches scipy.signal.windows.<name>(nperseg).
  k     <- seq(0L, nperseg - 1L)
  denom <- if (periodic) nperseg else nperseg - 1L
  win   <- switch(window,
    hann    = 0.5 * (1 - cos(2 * pi * k / denom)),
    hamming = 0.54 - 0.46 * cos(2 * pi * k / denom)
  )

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

  # Average across segments (vapply returns an n_one x n_segs matrix).
  # scipy's median is divided by .median_bias(n_segs) so that it estimates the
  # mean of the PSD for Gaussian noise.
  pgrams <- as.matrix(pgrams)
  psd <- if (average == "median") {
    rowmedian_cpp(pgrams) / .median_bias(ncol(pgrams))
  } else {
    rowMeans(pgrams)
  }

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
#' The defaults follow `yasa.bandpower()` (YASA 0.8.1): the PSD is
#' `scipy.signal.welch()` with a 4 s window (`win_sec = 4`), a periodic Hamming
#' window, median averaging over the overlapping segments (scipy's bias-corrected
#' median), an overlap of half a window, an FFT length equal to the window
#' length, and constant detrending. Band power is the Simpson integral of the
#' PSD, and the default bands are YASA's, including gamma at 30-40 Hz.
#' Relative power matches YASA's `relative = True`: each band is divided by the
#' integral of the PSD between the lowest and the highest band edge (0.5-40 Hz
#' by default). With bands that tile that range and edges on FFT bins this equals
#' the sum of the band powers; otherwise it is larger.
#'
#' To use another estimate, change `window` (`"hann"`), `average` (`"mean"`),
#' `noverlap` or `nfft`.
#'
#' @param psg An `mrpheus_psg` object from [mrpheus::prepare_psg()].
#' @param channels Character vector. Channel labels to include. If `NULL`
#'   (default), all non-bad EEG channels are used.
#' @param bands Named list of length-2 numeric vectors defining frequency bands
#'   in Hz. Default matches YASA's `bandpower()`:
#'   ```
#'   list(delta = c(0.5, 4), theta = c(4, 8), alpha = c(8, 12),
#'        sigma = c(12, 16), beta = c(16, 30), gamma = c(30, 40))
#'   ```
#' @param relative Logical. If `TRUE`, each band power is divided by the PSD
#'   integral between the lowest and highest band edge (dimensionless). Matches
#'   YASA `relative = True`. Default `FALSE`; note that YASA's default is `TRUE`.
#' @param win_sec Numeric. Welch window length in seconds. Default `4`.
#' @param noverlap Integer or `NULL`. Number of overlapping samples between
#'   Welch windows. `NULL` (default) uses half the window, as `scipy.signal.welch`
#'   does. Must be less than `win_sec * sample_rate`.
#' @param nfft Integer or `NULL`. FFT length. `NULL` (default) uses the window
#'   length in samples. If smaller than the window length it is automatically
#'   raised to the next power of 2.
#' @param window Character. Welch window, `"hamming"` (default, YASA's) or
#'   `"hann"`. Both are the periodic form that `scipy.signal.welch` uses.
#' @param average Character. How the segment periodograms are combined:
#'   `"median"` (default, YASA's; bias-corrected as in scipy) or `"mean"`.
#'
#' @return A tibble with columns `epoch`, `channel`, one column per named band,
#'   and `total_power` (the PSD integral between the lowest and highest band
#'   edge, in V^2 units before any relative scaling; YASA's `TotalAbsPow`).
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
#'
#' # The mean of Hann-windowed segments instead of YASA's defaults
#' bp <- compute_band_power(psg, window = "hann", average = "mean")
#' }
compute_band_power <- function(psg,
                               channels = NULL,
                               bands    = list(
                                 delta = c(0.5,  4),
                                 theta = c(4,    8),
                                 alpha = c(8,   12),
                                 sigma = c(12,  16),
                                 beta  = c(16,  30),
                                 gamma = c(30,  40)
                               ),
                               relative = FALSE,
                               win_sec  = 4,
                               noverlap = NULL,
                               nfft     = NULL,
                               window   = c("hamming", "hann"),
                               average  = c("median", "mean")) {
  stopifnot(inherits(psg, "mrpheus_psg"))
  window  <- match.arg(window)
  average <- match.arg(average)

  if (is.null(channels)) {
    channels <- psg$channel_map$label[
      psg$channel_map$type == "EEG" & !psg$channel_map$bad
    ]
  }

  # Frequency range over which the total power is integrated (YASA: the lowest
  # lower edge to the highest upper edge of the bands).
  fmin <- min(vapply(bands, `[`, numeric(1), 1L))
  fmax <- max(vapply(bands, `[`, numeric(1), 2L))

  rows <- lapply(seq_along(psg$epochs), function(i) {
    ep <- psg$epochs[[i]]

    lapply(channels, function(ch) {
      sig <- ep[[ch]]
      sr  <- psg$channel_map$sample_rate[psg$channel_map$label == ch]
      if (is.null(sig) || length(sig) < 2L) return(NULL)

      psd <- .welch_psd_bp(sig, fs = sr,
                            win_sec  = win_sec,
                            noverlap = noverlap,
                            nfft     = nfft,
                            window   = window,
                            average  = average,
                            periodic = TRUE)

      # As yasa.bandpower(): restrict the PSD to [fmin, fmax] first, then
      # integrate each band (and the whole range) with Simpson's rule.
      in_range <- psd$freq >= fmin & psd$freq <= fmax
      freq     <- psd$freq[in_range]
      spec     <- psd$spec[in_range]
      res      <- if (length(freq) >= 2L) freq[2L] - freq[1L] else NA_real_

      total <- if (is.na(res)) NA_real_ else .scipy_simpson(spec, dx = res)

      band_power <- vapply(bands, function(b) {
        idx <- freq >= b[1] & freq <= b[2]
        if (!any(idx) || is.na(res)) return(NA_real_)
        .scipy_simpson(spec[idx], dx = res)
      }, numeric(1))

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
