#' Detect slow oscillations
#'
#' Detects slow oscillations (SOs) in EEG channels using a zero-crossing
#' approach in the delta band (0.5-2 Hz), following the algorithm described in
#' Molle et al. (2002) and implemented in YASA (Vallat & Walker, 2021).
#'
#' @param psg An `mrpheus_psg` object from [mrpheus::prepare_psg()].
#' @param channel Character. EEG channel label. If `NULL` (default), the first
#'   non-bad EEG channel is used.
#' @param stages Integer vector or `NULL`. Epoch indices restricted to N2/N3.
#'   If `NULL`, detection runs across all epochs.
#' @param freq_so Numeric vector of length 2. SO frequency band (Hz).
#'   Default `c(0.5, 2)`.
#' @param amp_ptp_threshold_uv Numeric vector of length 2. Min and max
#'   acceptable peak-to-peak amplitude (uV). Default `c(75, 500)`.
#' @param duration_pos_s Numeric vector of length 2. Min and max positive half-
#'   wave duration (s). Default `c(0.1, 1.0)`.
#' @param duration_neg_s Numeric vector of length 2. Min and max negative half-
#'   wave duration (s). Default `c(0.1, 1.5)`.
#'
#' @return A tibble with one row per detected slow oscillation:
#' \describe{
#'   \item{epoch}{Integer.}
#'   \item{start_s}{Numeric. Onset (s) relative to epoch start.}
#'   \item{end_s}{Numeric. Offset (s) relative to epoch start.}
#'   \item{duration_s}{Numeric.}
#'   \item{neg_peak_uv}{Numeric. Negative peak amplitude (uV).}
#'   \item{pos_peak_uv}{Numeric. Positive peak amplitude (uV).}
#'   \item{ptp_uv}{Numeric. Peak-to-peak amplitude (uV).}
#'   \item{channel}{Character.}
#' }
#'
#' @references
#' Molle, M., Marshall, L., Gais, S., & Born, J. (2002). Grouping of spindle
#' activity during slow oscillations in human non-rapid eye movement sleep.
#' *Journal of Neuroscience*, 22(24), 10941-10947.
#'
#' Vallat, R., & Walker, M. P. (2021). An open-source, high-performance tool
#' for automated sleep staging. *eLife*, 10, e70092.
#' \doi{10.7554/eLife.70092}
#'
#' @export
compute_slow_oscillations <- function(psg,
                                      channel              = NULL,
                                      stages               = NULL,
                                      freq_so              = c(0.5, 2.0),
                                      amp_ptp_threshold_uv = c(75, 500),
                                      duration_pos_s       = c(0.1, 1.0),
                                      duration_neg_s       = c(0.1, 1.5)) {
  stopifnot(inherits(psg, "mrpheus_psg"))

  if (is.null(channel)) {
    channel <- psg$channel_map$label[
      psg$channel_map$type == "EEG" & !psg$channel_map$bad
    ][1]
  }

  sr        <- psg$channel_map$sample_rate[psg$channel_map$label == channel]
  epoch_idx <- stages %||% seq_along(psg$epochs)

  # Band-pass as a high-pass / low-pass cascade. A single 8th-order 0.5-2 Hz
  # band-pass is numerically singular in transfer-function form at 256 Hz
  # (gsignal::filter_zi: rcond ~ 4e-18), so each edge is designed separately.
  nyq   <- sr / 2
  bf_hp <- gsignal::butter(4, freq_so[1] / nyq, type = "high")
  bf_lp <- gsignal::butter(4, freq_so[2] / nyq, type = "low")

  sos <- lapply(epoch_idx, function(i) {
    sig <- psg$epochs[[i]][[channel]]
    if (is.null(sig) || length(sig) < 2) return(NULL)

    sig_filt <- as.numeric(
      gsignal::filtfilt(bf_lp, as.numeric(gsignal::filtfilt(bf_hp, sig)))
    )

    cands <- detect_so_candidates_cpp(
      sig_filt, sr,
      duration_neg_s[1], duration_neg_s[2],
      duration_pos_s[1], duration_pos_s[2],
      amp_ptp_threshold_uv[1], amp_ptp_threshold_uv[2]
    )
    if (nrow(cands) == 0L) return(NULL)

    tibble::tibble(
      epoch       = i,
      start_s     = cands[, 1L] / sr,
      end_s       = cands[, 3L] / sr,
      duration_s  = (cands[, 3L] - cands[, 1L]) / sr,
      neg_peak_uv = cands[, 4L],
      pos_peak_uv = cands[, 5L],
      ptp_uv      = cands[, 6L],
      channel     = channel
    )
  })

  out <- dplyr::bind_rows(sos)
  cli::cli_alert_success(
    "Detected {nrow(out)} slow oscillations in channel {.val {channel}}."
  )
  out
}
