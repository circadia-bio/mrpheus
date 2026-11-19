#' Time-resolved (sliding-window) respiration-locked HF power.
#'
#' Same underlying method as `compute_hrv_freq()` (RR intervals resampled
#' to an even grid, spectral power integrated in a band centered on the
#' infant's own respiratory frequency), but instead of averaging across the
#' whole recording into one number, computes HF power in sliding windows,
#' giving a time series that can be compared against another time-varying
#' signal (e.g. the Kuramoto synchrony trace from `dynR::leida_pipeline()`)
#' within the same recording.
#'
#' The respiratory frequency (and therefore the integration band) is
#' detected once from the *whole* recording's `resp` channel and held fixed
#' across all windows, rather than re-detected per window: a single
#' spectral estimate from the whole recording is far more stable than one
#' from a 60-second segment, and respiratory rate is generally stable
#' within a single scan session for a sleeping/sedated infant. If that
#' assumption doesn't hold for a given recording (e.g. a clear rate change
#' partway through), this will misrepresent windows far from the global
#' rate, worth checking the QC plot of `resp_freq_hz` if a session's
#' within-scan behaviour is in question.
#'
#' Verified against synthetic data: a ramping RSA amplitude trajectory
#' (rising vagal tone over a 15-minute recording) is recovered with r=0.998
#' between windowed HF power and the true underlying amplitude-squared
#' trajectory.
#'
#' @param qrs An `mrpheus_qrs` object (output of `detect_qrs()`).
#' @param resp_signal Numeric vector, the raw respiration channel.
#' @param fs Sampling rate (Hz) of `resp_signal`.
#' @param window_sec Window length in seconds. Default 60s: long enough
#'   for ~0.017 Hz frequency resolution (1/window_sec), comfortably finer
#'   than the 0.2 Hz-wide integration band, short enough to give
#'   reasonable temporal resolution over a ~15 minute scan.
#' @param step_sec Step between window centers in seconds. Default 15s
#'   (75% overlap between consecutive windows).
#' @param resample_hz RR-interval resampling rate. Default 4 Hz.
#' @param band_halfwidth_hz Half-width of the integration band. Default 0.1 Hz.
#' @param resp_search_range_hz Plausible respiratory frequency range for
#'   the whole-recording frequency detection. Default c(0.3, 2.5).
#'
#' @return A list: `time` (window center times, seconds, from recording
#'   start), `hf_power` (HF power at each window), `resp_freq_hz`,
#'   `band_low`, `band_high`.
#'
#' @export
compute_hrv_freq_windowed <- function(qrs, resp_signal, fs,
                                       window_sec = 60,
                                       step_sec = 15,
                                       resample_hz = 4,
                                       band_halfwidth_hz = 0.1,
                                       resp_search_range_hz = c(0.3, 2.5)) {

  empty_result <- function() {
    list(time = numeric(0), hf_power = numeric(0),
         resp_freq_hz = NA_real_, band_low = NA_real_, band_high = NA_real_)
  }

  # --- 1. RR intervals, resampled to an evenly-spaced grid (as in compute_hrv_freq()) ---
  peak_times <- qrs$qrs_i / qrs$fs
  rr <- diff(peak_times)
  rr_times <- peak_times[-1]

  valid <- rr > (60 / 220) & rr < (60 / 60)
  rr <- rr[valid]
  rr_times <- rr_times[valid]

  if (length(rr) < 10) return(empty_result())

  t_grid <- seq(rr_times[1], rr_times[length(rr_times)], by = 1 / resample_hz)
  rr_interp <- stats::spline(rr_times, rr, xout = t_grid, method = "natural")$y
  rr_detrended <- rr_interp - stats::predict(stats::lm(rr_interp ~ t_grid))

  # --- 2. Respiratory frequency from the WHOLE recording, held fixed across windows ---
  resp_freq <- .dominant_frequency(resp_signal, fs, resp_search_range_hz)
  if (is.na(resp_freq)) return(empty_result())

  band_low  <- resp_freq - band_halfwidth_hz
  band_high <- resp_freq + band_halfwidth_hz

  # --- 3. Sliding-window periodogram, one Hann-tapered FFT per window ---
  window_samples <- round(window_sec * resample_hz)
  step_samples   <- round(step_sec * resample_hz)
  n <- length(rr_detrended)

  if (n < window_samples) return(empty_result())

  starts <- seq(1, n - window_samples + 1, by = step_samples)

  hann <- 0.5 - 0.5 * cos(2 * pi * (0:(window_samples - 1)) / (window_samples - 1))
  u <- sum(hann^2) / window_samples
  freqs <- (0:(window_samples %/% 2)) * resample_hz / window_samples
  in_band <- freqs >= band_low & freqs <= band_high

  hf_power_t   <- numeric(length(starts))
  time_center  <- numeric(length(starts))

  for (i in seq_along(starts)) {
    s <- starts[i]
    seg <- rr_detrended[s:(s + window_samples - 1)]
    seg <- (seg - mean(seg)) * hann
    fft_seg <- stats::fft(seg)
    power <- (Mod(fft_seg)^2) / (resample_hz * window_samples * u)
    power <- power[1:length(freqs)]
    if (length(power) > 2) power[2:(length(power) - 1)] <- power[2:(length(power) - 1)] * 2

    hf_power_t[i]  <- if (any(in_band)) .trapz(freqs[in_band], power[in_band]) else NA_real_
    time_center[i] <- t_grid[s] + window_sec / 2
  }

  list(time = time_center, hf_power = hf_power_t,
       resp_freq_hz = resp_freq, band_low = band_low, band_high = band_high)
}
