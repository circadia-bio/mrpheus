# R/signal_filters.R
#
# Low-level signal filtering functions for physiological recordings.
# Each function operates on a plain numeric vector and can be called
# independently of the mrpheus_psg pipeline.
#
# Used internally by preprocess_psg() but exported so users can apply
# them in any context (e.g. filtering a single EEG channel loaded from
# a custom source).

# ── DC removal ───────────────────────────────────────────────────────────────

#' Remove DC offset from a signal
#'
#' Subtracts the signal mean (per-channel DC offset removal). Should be applied
#' before filtering to prevent filter transients from a non-zero baseline.
#'
#' @param signal Numeric vector.
#'
#' @return Numeric vector with mean removed, same length as `signal`.
#'
#' @export
#'
#' @examples
#' x <- c(100.1, 100.5, 99.8, 100.3)
#' remove_dc(x)
remove_dc <- function(signal) {
  signal - mean(signal, na.rm = TRUE)
}

# ── Powerline detection ───────────────────────────────────────────────────────

#' Auto-detect powerline frequency
#'
#' Estimates whether dominant powerline interference is at 50 Hz (Europe/UK)
#' or 60 Hz (North America) by comparing PSD power at both frequencies using
#' Welch's method on the first 10 seconds of the signal.
#'
#' If the signal's Nyquist limit is below 60 Hz, returns `50L` immediately.
#'
#' @param signal Numeric vector. Should be an EEG or broadband channel with
#'   sufficient bandwidth.
#' @param sr Numeric. Sampling rate in Hz. Must be > 100.
#'
#' @return Integer. `50L` or `60L`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' detect_powerline(eeg_signal, sr = 256)
#' }
detect_powerline <- function(signal, sr) {
  nyquist <- sr / 2

  if (60 >= nyquist) return(50L)

  n_use <- min(length(signal), as.integer(10 * sr))
  seg   <- signal[seq_len(n_use)]

  nfft <- 2L ^ ceiling(log2(4 * sr))
  psd  <- gsignal::pwelch(seg, fs = sr, window = nfft,
                           overlap = 0.5, nfft = nfft)

  idx_50 <- which.min(abs(psd$freq - 50))
  idx_60 <- which.min(abs(psd$freq - 60))

  if (psd$spec[idx_50] > psd$spec[idx_60] * 1.5) 50L else 60L
}

# ── Notch filter ──────────────────────────────────────────────────────────────

#' Apply a notch (bandstop) filter to a signal
#'
#' Designs a 2nd-order Butterworth bandstop filter centred at `freq` Hz with
#' bandwidth `bandwidth_hz` and applies it zero-phase using
#' `gsignal::filtfilt()`. Because the filter has only ~9 coefficients, it is
#' safe to use on long recordings without memory issues.
#'
#' By default also applies the filter at all harmonics of `freq` below the
#' Nyquist limit, which is the standard approach for powerline noise removal.
#'
#' @param signal Numeric vector.
#' @param sr Numeric. Sampling rate in Hz.
#' @param freq Numeric. Centre frequency in Hz (e.g. `50` or `60`).
#' @param bandwidth_hz Numeric. Full bandwidth of the notch in Hz. Default `2`.
#' @param harmonics Logical. Also notch harmonics (2×, 3×, ...) below
#'   Nyquist − 5 Hz. Default `TRUE`.
#'
#' @return Filtered numeric vector, same length as `signal`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' clean <- notch_filter(eeg_signal, sr = 256, freq = 50)
#' clean <- notch_filter(eeg_signal, sr = 256, freq = 60, harmonics = FALSE)
#' }
notch_filter <- function(signal, sr, freq, bandwidth_hz = 2, harmonics = TRUE) {
  nyquist <- sr / 2

  freqs <- freq
  if (harmonics) {
    harm <- 2 * freq
    while (harm < nyquist - 5) {
      freqs <- c(freqs, harm)
      harm  <- harm + freq
    }
  }

  for (f in freqs) {
    if (f >= nyquist) break
    half_bw <- bandwidth_hz / 2
    w_low   <- max((f - half_bw) / nyquist, 1e-4)
    w_high  <- min((f + half_bw) / nyquist, 1 - 1e-4)
    filt    <- gsignal::butter(2, c(w_low, w_high), type = "stop")
    signal  <- as.numeric(gsignal::filtfilt(filt, signal))
  }

  signal
}

# ── Bandpass filter ───────────────────────────────────────────────────────────

#' Apply a bandpass filter to a signal
#'
#' Designs a 4th-order Butterworth bandpass filter and applies it zero-phase
#' using `gsignal::filtfilt()`. Handles edge cases gracefully:
#'
#' - `low_hz <= 0` → low-pass at `high_hz`
#' - `high_hz >= Nyquist` → high-pass at `low_hz`
#' - Otherwise → bandpass between `low_hz` and `high_hz`
#'
#' The `high_hz` cutoff is clamped to 99 % of Nyquist to avoid instability.
#' Because this uses a short IIR filter, it is safe on long recordings.
#'
#' @param signal Numeric vector.
#' @param sr Numeric. Sampling rate in Hz.
#' @param low_hz Numeric. Low cutoff frequency in Hz. Use `0` or negative to
#'   apply a low-pass filter only.
#' @param high_hz Numeric. High cutoff frequency in Hz.
#'
#' @return Filtered numeric vector, same length as `signal`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Standard EEG bandpass
#' filtered <- bandpass_filter(eeg_signal, sr = 256, low_hz = 0.3, high_hz = 35)
#'
#' # High-pass only (baseline drift removal)
#' hp <- bandpass_filter(eeg_signal, sr = 256, low_hz = 0.1, high_hz = Inf)
#' }
bandpass_filter <- function(signal, sr, low_hz, high_hz) {
  nyquist <- sr / 2
  high_hz <- min(high_hz, nyquist * 0.99)

  filt <- if (low_hz <= 0) {
    gsignal::butter(4, high_hz / nyquist, type = "low")
  } else if (high_hz >= nyquist * 0.99) {
    gsignal::butter(4, low_hz / nyquist, type = "high")
  } else {
    gsignal::butter(4, c(low_hz, high_hz) / nyquist, type = "pass")
  }

  as.numeric(gsignal::filtfilt(filt, signal))
}
