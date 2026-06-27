#' Export a staged hypnogram to hypnor
#'
#' Converts the staging tibble produced by [stage_epochs()] into a
#' `hypnor_hypnogram` object ready for architecture metric computation and
#' visualisation in the `hypnor` package. This is the primary handoff between
#' `mrpheus` and the rest of the Circadia Lab ecosystem.
#'
#' @param staging A tibble from [stage_epochs()] with columns `epoch`, `stage`,
#'   and optional probability columns.
#' @param epoch_s Numeric. Epoch duration in seconds. Must match the `epoch_s`
#'   used in [prepare_psg()]. Default `30`.
#' @param start_time POSIXct or `NULL`. Recording start time. Used to compute
#'   clock-time axes in `hypnor` plots. If `NULL`, epochs are indexed from 0.
#' @param participant_id Character or `NULL`. Optional identifier passed through
#'   to `hypnor` and `syncR`.
#'
#' @return An object of class `hypnor_hypnogram` (defined in `hypnor`).
#'   If `hypnor` is not installed, returns the staging tibble invisibly with a
#'   message.
#'
#' @seealso [stage_epochs()]
#'
#' @export
export_hypnogram <- function(staging,
                             epoch_s        = 30,
                             start_time     = NULL,
                             participant_id = NULL) {
  if (!is.data.frame(staging) || !"stage" %in% names(staging)) {
    cli::cli_abort("{.arg staging} must be a tibble from {.fn stage_epochs}.")
  }

  if (!requireNamespace("hypnor", quietly = TRUE)) {
    cli::cli_alert_warning(
      "{.pkg hypnor} is not installed. Returning staging tibble as-is."
    )
    return(invisible(staging))
  }

  hypnor::new_hypnogram(
    stages         = staging$stage,
    epoch_s        = epoch_s,
    start_time     = start_time,
    participant_id = participant_id,
    source         = "mrpheus",
    resolution     = "AASM"
  )
}
