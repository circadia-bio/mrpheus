#' Export a staged hypnogram for use with hypnor
#'
#' Prepares the staging tibble produced by [mrpheus::stage_epochs()] for
#' downstream use with the `hypnor` package. Attaches recording metadata as
#' attributes and returns a tibble of class `mrpheus_hypnogram`, which
#' `hypnor::new_hypnogram()` accepts directly once `hypnor` is installed.
#'
#' @param staging A tibble from [mrpheus::stage_epochs()] with columns `epoch`,
#'   `stage`, and optional probability columns.
#' @param epoch_s Numeric. Epoch duration in seconds. Must match the `epoch_s`
#'   used in [mrpheus::prepare_psg()]. Default `30`.
#' @param start_time POSIXct or `NULL`. Recording start time. Used to compute
#'   clock-time axes in `hypnor` plots. If `NULL`, epochs are indexed from 0.
#' @param participant_id Character or `NULL`. Optional identifier passed through
#'   to `hypnor` and `syncR`.
#'
#' @return A tibble of class `mrpheus_hypnogram` with columns `epoch`, `stage`,
#'   and any probability columns from the staging model. Metadata (`epoch_s`,
#'   `start_time`, `participant_id`, `source`, `resolution`) are attached as
#'   attributes and forwarded to `hypnor::new_hypnogram()` when `hypnor` is
#'   available.
#'
#' @seealso [mrpheus::stage_epochs()]
#'
#' @export
export_hypnogram <- function(staging,
                             epoch_s        = 30,
                             start_time     = NULL,
                             participant_id = NULL) {
  if (!is.data.frame(staging) || !"stage" %in% names(staging)) {
    cli::cli_abort("`staging` must be a tibble from `stage_epochs()`.")
  }

  out <- staging
  attr(out, "epoch_s")        <- epoch_s
  attr(out, "start_time")     <- start_time
  attr(out, "participant_id") <- participant_id
  attr(out, "source")         <- "mrpheus"
  attr(out, "resolution")     <- "AASM"
  class(out) <- c("mrpheus_hypnogram", class(out))

  cli::cli_alert_success(
    "Hypnogram ready: {nrow(out)} epochs. \\
    Pass to {.code hypnor::new_hypnogram()} once {.pkg hypnor} is available."
  )
  out
}

#' @export
print.mrpheus_hypnogram <- function(x, ...) {
  cli::cli_h1("mrpheus hypnogram")
  cli::cli_inform(c(
    "i" = "Epochs: {nrow(x)}",
    "i" = "Epoch length: {attr(x, 'epoch_s')} s",
    "i" = "Participant: {attr(x, 'participant_id') %||% 'unset'}",
    "i" = "Source: {attr(x, 'source')} / {attr(x, 'resolution')}"
  ))
  NextMethod()
  invisible(x)
}
