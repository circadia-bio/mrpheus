#' Detect sleep spindles
#'
#' Detects sleep spindles in EEG channels using a band-pass / RMS envelope
#' approach, closely following the algorithm in Lacourse et al. (2019) and
#' implemented in YASA (Vallat & Walker, 2021). Spindles are identified as
#' transient bursts of 11–16 Hz activity during NREM sleep.
#'
#' @param psg An `mrpheus_psg` object from [mrpheus::prepare_psg()].
#' @param channel Character. EEG channel label. If `NULL` (default), the first
#'   non-bad EEG channel is used.
#' @param stages Integer vector or `NULL`. Epoch indices (1-based) restricted to
#'   NREM sleep. If `NULL`, spindles are detected across all epochs.
#' @param freq_spindle Numeric vector of length 2. Spindle frequency band (Hz).
#'   Default `c(11, 16)`.
#' @param rms_window_s Numeric. Moving RMS window length in seconds. Default `0.3`.
#' @param min_duration_s Numeric. Minimum spindle duration in seconds.
#'   Default `0.5`.
#' @param max_duration_s Numeric. Maximum spindle duration in seconds.
#'   Default `3.0`.
#' @param threshold_sd Numeric. RMS threshold as a multiple of the channel SD
#'   (within sigma band). Default `1.5`.
#'
#' @return A tibble with one row per detected spindle:
#' \describe{
#'   \item{epoch}{Integer. Epoch in which the spindle was detected.}
#'   \item{start_s}{Numeric. Spindle onset relative to epoch start (seconds).}
#'   \item{end_s}{Numeric. Spindle offset relative to epoch start (seconds).}
#'   \item{duration_s}{Numeric. Spindle duration in seconds.}
#'   \item{peak_freq_hz}{Numeric. Peak frequency within the spindle.}
#'   \item{rms_uv}{Numeric. Mean RMS amplitude within the spindle (µV).}
#'   \item{channel}{Character. Channel label.}
#' }
#'
#' @references
#' Lacourse, K., Yetton, B., Mednick, S., & Bhatt, D. L. (2019).
#' *Massive online sleep staging: A polysomnography data repository*.
#' Journal of Sleep Research.
#'
#' Vallat, R., & Walker, M. P. (2021). An open-source, high-performance tool
#' for automated sleep staging. *eLife*, 10, e70092.
#' \doi{10.7554/eLife.70092}
#'
#' @export
compute_spindles <- function(psg,
                             channel        = NULL,
                             stages         = NULL,
                             freq_spindle   = c(11, 16),
                             rms_window_s   = 0.3,
                             min_duration_s = 0.5,
                             max_duration_s = 3.0,
                             threshold_sd   = 1.5) {
  stopifnot(inherits(psg, "mrpheus_psg"))

  if (is.null(channel)) {
    channel <- psg$channel_map$label[
      psg$channel_map$type == "EEG" & !psg$channel_map$bad
    ][1]
  }

  sr        <- psg$channel_map$sample_rate[psg$channel_map$label == channel]
  epoch_idx <- stages %||% seq_along(psg$epochs)

  spindles <- lapply(epoch_idx, function(i) {
    sig <- psg$epochs[[i]][[channel]]
    if (is.null(sig) || length(sig) < 2) return(NULL)

    bf  <- gsignal::butter(4, freq_spindle / (sr / 2), type = "pass")
    sig_filt <- gsignal::filtfilt(bf, sig)

    win_samples <- round(rms_window_s * sr)
    rms_env <- zoo::rollapply(sig_filt^2, width = win_samples,
                              FUN = function(x) sqrt(mean(x)),
                              fill = NA, align = "center")

    threshold <- threshold_sd * sd(sig_filt, na.rm = TRUE)
    above     <- !is.na(rms_env) & rms_env > threshold

    if (!any(above)) return(NULL)

    runs  <- rle(above)
    ends  <- cumsum(runs$lengths)
    starts <- ends - runs$lengths + 1L

    sp_list <- lapply(which(runs$values), function(k) {
      s <- starts[k]
      e <- ends[k]
      dur_s <- (e - s) / sr
      if (dur_s < min_duration_s || dur_s > max_duration_s) return(NULL)

      seg <- sig_filt[s:e]
      if (length(seg) < 4) return(NULL)
      psd      <- gsignal::pwelch(seg, fs = sr,
                                   nfft = max(64L, 2^ceiling(log2(length(seg)))),
                                   window = length(seg))
      peak_idx <- which.max(psd$spec)

      tibble::tibble(
        epoch        = i,
        start_s      = s / sr,
        end_s        = e / sr,
        duration_s   = dur_s,
        peak_freq_hz = psd$freq[peak_idx],
        rms_uv       = mean(rms_env[s:e], na.rm = TRUE),
        channel      = channel
      )
    })
    dplyr::bind_rows(sp_list)
  })

  out <- dplyr::bind_rows(spindles)
  cli::cli_alert_success("Detected {nrow(out)} spindles in channel {.val {channel}}.")
  out
}
