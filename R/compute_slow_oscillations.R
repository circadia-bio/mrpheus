#' Detect slow oscillations
#'
#' Detects slow oscillations (SOs) in EEG channels using a zero-crossing
#' approach, following the algorithm described in Molle et al. (2002). The
#' defaults (band, duration, trough, peak and peak-to-peak limits) are those of
#' `yasa.sw_detect()` in YASA (Vallat & Walker, 2021), so that densities are
#' comparable in magnitude. The detection itself is not identical to YASA's:
#' this function pairs zero crossings into a negative and a following positive
#' half-wave, whereas YASA pairs each negative peak with the nearest following
#' positive peak, so event counts and timings differ. Epochs are filtered one
#' at a time.
#'
#' @param psg An `mrpheus_psg` object from [mrpheus::prepare_psg()].
#' @param channel Character. EEG channel label. If `NULL` (default), the first
#'   non-bad EEG channel is used.
#' @param stages Integer vector or `NULL`. Epoch indices restricted to N2/N3.
#'   If `NULL`, detection runs across all epochs.
#' @param freq_so Numeric vector of length 2. SO frequency band (Hz), as in
#'   YASA's `freq_sw`. Default `c(0.3, 1.5)`. See `transition_hz`.
#' @param amp_ptp_threshold_uv Numeric vector of length 2. Min and max
#'   acceptable peak-to-peak amplitude (uV). Default `c(75, 350)`.
#' @param duration_pos_s Numeric vector of length 2. Min and max positive half-
#'   wave duration (s). Default `c(0.1, 1.0)`.
#' @param duration_neg_s Numeric vector of length 2. Min and max negative half-
#'   wave duration (s). Default `c(0.3, 1.5)`.
#' @param amp_neg_uv Numeric vector of length 2, or `NULL`. Min and max absolute
#'   amplitude of the negative peak (uV). Default `c(40, 200)`. `NULL` applies
#'   no limit.
#' @param amp_pos_uv Numeric vector of length 2, or `NULL`. Min and max
#'   amplitude of the positive peak (uV). Default `c(10, 150)`. `NULL` applies
#'   no limit.
#' @param transition_hz Numeric. Width of the filter transition band (Hz).
#'   YASA band-passes with an FIR filter whose -6 dB points lie half a
#'   transition band outside `freq_so` (0.2 and 1.6 Hz for the defaults). Here
#'   the band-pass is a zero-phase Butterworth high-pass / low-pass cascade
#'   whose -6 dB points are placed at `freq_so[1] - transition_hz / 2` and
#'   `freq_so[2] + transition_hz / 2`. Default `0.2`.
#'
#' @details
#' Amplitudes are assumed to be in uV. The high-pass is a 4th-order Butterworth
#' in transfer-function form; its coefficients become ill-conditioned at very
#' low cutoffs and high sample rates (rcond about 2e-14 at 1000 Hz for a 0.2 Hz
#' cutoff, against a failure limit near 2e-16), so very high sample rates
#' should be downsampled first.
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
                                      freq_so              = c(0.3, 1.5),
                                      amp_ptp_threshold_uv = c(75, 350),
                                      duration_pos_s       = c(0.1, 1.0),
                                      duration_neg_s       = c(0.3, 1.5),
                                      amp_neg_uv           = c(40, 200),
                                      amp_pos_uv           = c(10, 150),
                                      transition_hz        = 0.2) {
  stopifnot(inherits(psg, "mrpheus_psg"))

  if (is.null(channel)) {
    channel <- psg$channel_map$label[
      psg$channel_map$type == "EEG" & !psg$channel_map$bad
    ][1]
  }

  sr        <- psg$channel_map$sample_rate[psg$channel_map$label == channel]
  epoch_idx <- stages %||% seq_along(psg$epochs)

  # Band-pass as a high-pass / low-pass cascade. A single 8th-order narrow
  # band-pass is numerically singular in transfer-function form at 256 Hz
  # (gsignal::filter_zi: rcond ~ 4e-18 for 0.5-2 Hz), so each edge is designed
  # separately. The -6 dB points (the cutoffs, after the zero-phase double
  # pass) sit half a transition band outside freq_so, as in YASA's FIR filter.
  nyq    <- sr / 2
  hp_cut <- freq_so[1] - transition_hz / 2
  lp_cut <- freq_so[2] + transition_hz / 2
  if (hp_cut <= 0 || lp_cut >= nyq || hp_cut >= lp_cut) {
    cli::cli_abort(c(
      "Invalid slow oscillation filter band.",
      "i" = "High-pass {hp_cut} Hz and low-pass {lp_cut} Hz must satisfy 0 < high-pass < low-pass < {nyq} Hz.",
      "i" = "These are `freq_so` widened by `transition_hz / 2` on each side."
    ))
  }
  bf_hp <- gsignal::butter(4, hp_cut / nyq, type = "high")
  bf_lp <- gsignal::butter(4, lp_cut / nyq, type = "low")

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

    # Trough and peak amplitude limits (YASA's amp_neg / amp_pos). They do not
    # change how crossings are paired, so applying them here is equivalent to
    # applying them inside the C++ walk.
    if (nrow(cands) > 0L && !is.null(amp_neg_uv))
      cands <- cands[-cands[, 4L] >= amp_neg_uv[1] & -cands[, 4L] <= amp_neg_uv[2],
                     , drop = FALSE]
    if (nrow(cands) > 0L && !is.null(amp_pos_uv))
      cands <- cands[cands[, 5L] >= amp_pos_uv[1] & cands[, 5L] <= amp_pos_uv[2],
                     , drop = FALSE]
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
