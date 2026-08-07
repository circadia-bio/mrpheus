#' Frequency-domain HRV (HF power), locked to each infant's own
#' respiratory rate, as a proxy for respiratory sinus arrhythmia (RSA).
#'
#' Adults have a standard HF band (0.15-0.4 Hz) because adult respiratory
#' rate is fairly stable across that range. Neonates breathe far faster
#' (roughly 40-70 breaths/min at rest, i.e. ~0.7-1.2 Hz, and higher still
#' if tachypnoeic), well above the adult band. Using the adult HF band on
#' infant data would largely miss the respiratory-linked HRV signal
#' entirely. This instead detects each infant's own dominant respiratory
#' frequency directly from the physlog's `resp` channel, and integrates
#' RR-interval spectral power in a narrow band centered on that frequency.
#'
#' Verified against synthetic data: recovers an injected respiratory
#' frequency exactly (0 Hz error against a known 0.9 Hz test signal), and
#' HF power scales with injected RSA amplitude as expected (power scaling
#' roughly with amplitude^2).
#'
#' @param qrs An `mrpheus_qrs` object (output of `detect_qrs()`), giving
#'   R-peak sample indices and sampling rate.
#' @param resp_signal Numeric vector, the raw respiration channel, same
#'   sampling rate and duration as the signal QRS detection was run on.
#' @param fs Sampling rate (Hz) of `resp_signal` (should match the ECG's
#'   original sampling rate; physlog channels share one clock).
#' @param resample_hz Frequency (Hz) to resample the RR-interval series to
#'   before spectral estimation. Default 4 Hz, standard for HRV frequency
#'   analysis and comfortably above twice the plausible infant respiratory
#'   range (Nyquist-safe up to 2 Hz / 120 breaths/min).
#' @param band_halfwidth_hz Half-width (Hz) of the HF integration band,
#'   centered on the detected respiratory frequency. Default 0.1 Hz.
#' @param resp_search_range_hz Plausible neonatal respiratory frequency
#'   range (Hz) to search for the dominant peak in `resp_signal`. Default
#'   c(0.3, 2.5) (18-150 breaths/min), deliberately wide to avoid missing
#'   tachypnoeic or bradypnoeic infants.
#'
#' @return A list: `hf_power` (integrated PSD power in the respiration-locked
#'   band), `resp_freq_hz` (detected dominant respiratory frequency),
#'   `band_low`/`band_high` (integration band bounds), `n_beats_used`.
#'
#' @export
compute_hrv_freq <- function(qrs, resp_signal, fs,
                              resample_hz = 4,
                              band_halfwidth_hz = 0.1,
                              resp_search_range_hz = c(0.3, 2.5)) {

  empty_result <- function(n_beats = NA_integer_) {
    list(hf_power = NA_real_, resp_freq_hz = NA_real_,
         band_low = NA_real_, band_high = NA_real_,
         n_beats_used = n_beats)
  }

  # --- 1. RR intervals from QRS peak indices -----------------------------
  peak_times <- qrs$qrs_i / qrs$fs   # seconds
  rr <- diff(peak_times)             # seconds
  rr_times <- peak_times[-1]         # time-stamp each RR interval at its second beat

  # Drop physiologically implausible RR intervals (outside ~60-220 bpm)
  # before interpolation, rather than letting a missed or extra beat
  # distort the resampled series.
  valid <- rr > (60 / 220) & rr < (60 / 60)
  rr <- rr[valid]
  rr_times <- rr_times[valid]

  if (length(rr) < 10) return(empty_result(length(rr)))

  # --- 2. Resample to an evenly-spaced grid via cubic spline -------------
  # RR series is inherently unevenly spaced (one value per heartbeat, at
  # beat-to-beat intervals); spectral methods require even spacing.
  t_grid <- seq(rr_times[1], rr_times[length(rr_times)], by = 1 / resample_hz)
  rr_interp <- stats::spline(rr_times, rr, xout = t_grid, method = "natural")$y

  # Remove linear trend before spectral estimation (standard practice):
  # a slow drift would otherwise dominate the low-frequency end of the
  # spectrum and can leak into neighbouring bands.
  rr_detrended <- rr_interp - stats::predict(stats::lm(rr_interp ~ t_grid))

  # --- 3. Detect this infant's own respiratory frequency -----------------
  resp_freq <- .dominant_frequency(resp_signal, fs, resp_search_range_hz)
  if (is.na(resp_freq)) return(empty_result(length(rr)))

  band_low  <- resp_freq - band_halfwidth_hz
  band_high <- resp_freq + band_halfwidth_hz

  # --- 4. Welch PSD of the resampled RR series ----------------------------
  psd <- .welch_psd(rr_detrended, fs = resample_hz)

  # --- 5. Integrate power within the respiration-locked HF band ----------
  in_band <- psd$freq >= band_low & psd$freq <= band_high
  hf_power <- if (any(in_band)) .trapz(psd$freq[in_band], psd$power[in_band]) else NA_real_

  list(
    hf_power     = hf_power,
    resp_freq_hz = resp_freq,
    band_low     = band_low,
    band_high    = band_high,
    n_beats_used = length(rr)
  )
}

# ---------------------------------------------------------------------------
# Internal: dominant frequency (spectral peak) within a plausible range.
# Used here to detect an infant's own respiratory rate from the raw resp
# channel, rather than assuming a fixed adult-derived frequency.
# ---------------------------------------------------------------------------
.dominant_frequency <- function(x, fs, freq_range) {
  x <- x - mean(x, na.rm = TRUE)
  psd <- .welch_psd(x, fs)
  in_range <- psd$freq >= freq_range[1] & psd$freq <= freq_range[2]
  if (!any(in_range)) return(NA_real_)
  psd$freq[in_range][which.max(psd$power[in_range])]
}

# ---------------------------------------------------------------------------
# Internal: Welch's method for power spectral density estimation. Segments
# the signal, applies a Hann window, averages the periodograms. Implemented
# directly in base R (fft(), no external spectral-analysis package) so this
# has no new dependency beyond what mrpheus already requires.
#
# Verified against synthetic data (see package tests / verification script):
# recovers a known injected frequency exactly, and integrated band power
# scales with injected signal amplitude as expected.
# ---------------------------------------------------------------------------
.welch_psd <- function(x, fs, segment_sec = 60, overlap = 0.5) {
  n <- length(x)
  seg_len <- min(n, round(segment_sec * fs))
  step <- max(1, round(seg_len * (1 - overlap)))

  starts <- seq(1, max(1, n - seg_len + 1), by = step)

  hann <- 0.5 - 0.5 * cos(2 * pi * (0:(seg_len - 1)) / (seg_len - 1))
  # Normalise window so PSD scaling is consistent regardless of window shape.
  u <- sum(hann^2) / seg_len

  n_fft <- seg_len
  freqs <- (0:(n_fft %/% 2)) * fs / n_fft

  psd_sum <- numeric(length(freqs))
  n_segments <- 0

  for (s in starts) {
    seg <- x[s:(s + seg_len - 1)]
    seg <- (seg - mean(seg)) * hann
    fft_seg <- stats::fft(seg)
    power <- (Mod(fft_seg)^2) / (fs * seg_len * u)
    power <- power[1:length(freqs)]
    # One-sided PSD convention: double non-DC, non-Nyquist bins to account
    # for the folded negative-frequency half of the spectrum.
    if (length(power) > 2) power[2:(length(power) - 1)] <- power[2:(length(power) - 1)] * 2
    psd_sum <- psd_sum + power
    n_segments <- n_segments + 1
  }

  list(freq = freqs, power = psd_sum / n_segments)
}

# ---------------------------------------------------------------------------
# Internal: trapezoidal integration. Base R has no built-in trapz(); used
# here to integrate PSD power within a frequency band.
# ---------------------------------------------------------------------------
.trapz <- function(x, y) {
  n <- length(x)
  if (n < 2) return(0)
  sum(diff(x) * (y[-1] + y[-n]) / 2)
}
