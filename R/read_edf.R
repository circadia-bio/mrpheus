#' Read an EDF or EDF+ recording
#'
#' Reads a European Data Format (EDF or EDF+) file and returns a structured
#' list containing signal data, channel metadata, and recording header.
#' Wraps `edfReader::readEdfHeader()` and `edfReader::readEdfSignals()` with
#' consistent output formatting for the mrpheus pipeline.
#'
#' @param path Character. Path to an `.edf` or `.edf+` file.
#' @param channels Character vector or `NULL`. Channel labels to import.
#'   If `NULL` (default), all channels are imported.
#' @param only_header Logical. If `TRUE`, return only the header without
#'   reading signal data. Useful for quick channel inspection. Default `FALSE`.
#'
#' @return A list of class `mrpheus_edf` with components:
#' \describe{
#'   \item{header}{Data frame. Recording metadata (patient info, start time,
#'     number of signals, etc.).}
#'   \item{signals}{Named list of numeric vectors, one per channel.}
#'   \item{channels}{Data frame. Channel-level metadata: label, sample rate,
#'     physical min/max, digital min/max, transducer type, prefiltering.}
#'   \item{duration_s}{Numeric. Total recording duration in seconds.}
#'   \item{path}{Character. Resolved path to the source file.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' rec <- read_edf("data/psg_001.edf")
#' rec <- read_edf("data/psg_001.edf", channels = c("EEG Fpz-Cz", "EOG horizontal"))
#' rec$channels
#' }
read_edf <- function(path, channels = NULL, only_header = FALSE) {
  if (!file.exists(path)) {
    cli::cli_abort("File not found: {.path {path}}")
  }

  hdr <- edfReader::readEdfHeader(path)

  if (only_header) {
    return(structure(list(header = hdr, path = normalizePath(path)),
                     class = "mrpheus_edf_header"))
  }

  sigs <- edfReader::readEdfSignals(hdr, signals = channels %||% "All")

  # Build channel metadata table
  ch_meta <- tibble::tibble(
    label       = hdr$sHeaders$label,
    sample_rate = hdr$sHeaders$sRate,
    phys_min    = hdr$sHeaders$physicalMinimum,
    phys_max    = hdr$sHeaders$physicalMaximum,
    transducer  = hdr$sHeaders$transducerType,
    prefilter   = hdr$sHeaders$prefiltering
  )

  if (!is.null(channels)) {
    ch_meta <- ch_meta[ch_meta$label %in% channels, ]
  }

  structure(
    list(
      header     = hdr,
      signals    = sigs,
      channels   = ch_meta,
      duration_s = as.numeric(hdr$endTime - hdr$startTime),
      path       = normalizePath(path)
    ),
    class = "mrpheus_edf"
  )
}

#' @export
print.mrpheus_edf <- function(x, ...) {
  cli::cli_h1("mrpheus EDF recording")
  cli::cli_inform(c(
    "i" = "Path: {.path {x$path}}",
    "i" = "Duration: {round(x$duration_s / 3600, 2)} hours",
    "i" = "Channels: {nrow(x$channels)}"
  ))
  print(x$channels[, c("label", "sample_rate", "transducer")], ...)
  invisible(x)
}
