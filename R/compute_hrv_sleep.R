#' Compute HRV metrics across sleep stages
#'
#' Extracts R-peak positions from the ECG channel, computes standard time- and
#' frequency-domain heart rate variability (HRV) metrics, and stratifies results
#' by sleep stage.
#'
#' @param psg An `mrpheus_psg` object from [mrpheus::prepare_psg()].
#' @param staging Tibble from [mrpheus::stage_epochs()] or `NULL`. If `NULL`,
#'   HRV is computed across all epochs without stage stratification.
#' @param ecg_channel Character or `NULL`. ECG channel label. If `NULL`
#'   (default), the first non-bad ECG channel is used.
#' @param min_rr_ms Numeric. Minimum physiologically plausible RR interval (ms).
#'   Default `300` (200 bpm ceiling).
#' @param max_rr_ms Numeric. Maximum physiologically plausible RR interval (ms).
#'   Default `2000` (30 bpm floor).
#'
#' @return A tibble with one row per sleep stage (or one row if `staging` is
#'   `NULL`):
#' \describe{
#'   \item{stage}{Character. AASM stage or `"ALL"`.}
#'   \item{n_epochs}{Integer. Number of epochs in this stage.}
#'   \item{mean_rr_ms}{Numeric. Mean RR interval (ms).}
#'   \item{sdnn_ms}{Numeric. SDNN — SD of all NN intervals.}
#'   \item{rmssd_ms}{Numeric. RMSSD — root mean square of successive differences.}
#'   \item{lf_power}{Numeric. LF band power (0.04–0.15 Hz).}
#'   \item{hf_power}{Numeric. HF band power (0.15–0.4 Hz).}
#'   \item{lf_hf_ratio}{Numeric. LF/HF ratio.}
#' }
#'
#' @export
compute_hrv_sleep <- function(psg,
                              staging     = NULL,
                              ecg_channel = NULL,
                              min_rr_ms   = 300,
                              max_rr_ms   = 2000) {
  stopifnot(inherits(psg, "mrpheus_psg"))

  ecg_channel <- ecg_channel %||%
    psg$channel_map$label[psg$channel_map$type == "ECG" & !psg$channel_map$bad][1]

  if (is.na(ecg_channel)) {
    cli::cli_abort("No ECG channel available. Check `prepare_psg()` channel map.")
  }

  cli::cli_alert_warning(
    "`compute_hrv_sleep()` is a stub. Full R-peak detection implementation pending."
  )

  tibble::tibble(
    stage       = "ALL",
    n_epochs    = psg$n_epochs,
    mean_rr_ms  = NA_real_,
    sdnn_ms     = NA_real_,
    rmssd_ms    = NA_real_,
    lf_power    = NA_real_,
    hf_power    = NA_real_,
    lf_hf_ratio = NA_real_
  )
}
