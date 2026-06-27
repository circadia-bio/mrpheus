#' Detect respiratory events (apneas and hypopneas)
#'
#' Detects apneas and hypopneas in respiratory airflow and effort signals,
#' returning event-level metadata and summary indices. Detection follows AASM
#' 2012/2017 criteria: >= 90% signal reduction for >= 10 s (apnea), or >= 30%
#' reduction with >= 3% SpO2 desaturation or arousal for >= 10 s (hypopnea).
#'
#' @param psg An `mrpheus_psg` object from [mrpheus::prepare_psg()].
#' @param airflow_channel Character. Airflow channel label (e.g. nasal pressure
#'   or thermistor). If `NULL` (default), the first RESP channel is used.
#' @param spo2_channel Character or `NULL`. SpO2 channel for hypopnea scoring.
#' @param min_duration_s Numeric. Minimum event duration in seconds. Default `10`.
#' @param apnea_threshold Numeric. Fractional reduction required for apnea
#'   classification. Default `0.90`.
#' @param hypopnea_threshold Numeric. Fractional reduction required for
#'   hypopnea classification. Default `0.30`.
#' @param desaturation_threshold Numeric. SpO2 drop (percentage points) required
#'   to confirm hypopnea. Default `3`.
#'
#' @return A list with:
#' \describe{
#'   \item{events}{Tibble. One row per event: `epoch`, `start_s`, `end_s`,
#'     `duration_s`, `type` (`"apnea"` / `"hypopnea"`), `desaturation`.}
#'   \item{summary}{Tibble. `n_apneas`, `n_hypopneas`, `n_events`.}
#' }
#'
#' @seealso [mrpheus::compute_ahi()], [mrpheus::compute_odi()]
#'
#' @export
detect_apneas <- function(psg,
                          airflow_channel        = NULL,
                          spo2_channel           = NULL,
                          min_duration_s         = 10,
                          apnea_threshold        = 0.90,
                          hypopnea_threshold     = 0.30,
                          desaturation_threshold = 3) {
  stopifnot(inherits(psg, "mrpheus_psg"))

  airflow_channel <- airflow_channel %||%
    psg$channel_map$label[psg$channel_map$type == "RESP" & !psg$channel_map$bad][1]

  if (is.na(airflow_channel)) {
    cli::cli_abort("No RESP channel available. Check `prepare_psg()` channel map.")
  }

  cli::cli_alert_warning(
    "`detect_apneas()` is a stub. Full implementation pending."
  )

  list(
    events  = tibble::tibble(
      epoch = integer(), start_s = numeric(), end_s = numeric(),
      duration_s = numeric(), type = character(), desaturation = numeric()
    ),
    summary = tibble::tibble(n_apneas = 0L, n_hypopneas = 0L, n_events = 0L)
  )
}

#' Compute Apnea-Hypopnea Index (AHI)
#'
#' Calculates the AHI from respiratory event data and total sleep time.
#'
#' @param events Output of [mrpheus::detect_apneas()].
#' @param tst_hours Numeric. Total sleep time in hours.
#'
#' @return Numeric scalar. AHI (events per hour).
#'
#' @export
compute_ahi <- function(events, tst_hours) {
  stopifnot(is.list(events), "events" %in% names(events))
  nrow(events$events) / tst_hours
}

#' Compute Oxygen Desaturation Index (ODI)
#'
#' Calculates the ODI (number of >= 3% SpO2 desaturations per hour of sleep)
#' from the SpO2 channel.
#'
#' @param psg An `mrpheus_psg` object from [mrpheus::prepare_psg()].
#' @param spo2_channel Character. SpO2 channel label.
#' @param tst_hours Numeric. Total sleep time in hours.
#' @param threshold Numeric. Desaturation threshold (percentage points).
#'   Default `3`.
#'
#' @return Numeric scalar. ODI (desaturations per hour).
#'
#' @export
compute_odi <- function(psg, spo2_channel, tst_hours, threshold = 3) {
  stopifnot(inherits(psg, "mrpheus_psg"))
  cli::cli_alert_warning("`compute_odi()` is a stub. Full implementation pending.")
  NA_real_
}
