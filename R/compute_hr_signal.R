#' Compute Instantaneous Heart Rate from R-Peak Indices
#'
#' Converts R-peak sample indices (from [mrpheus::detect_qrs()]) into a
#' sample-by-sample instantaneous heart rate trace by filling each inter-beat
#' interval with the constant HR derived from that interval's length.
#'
#' @param qrs A `mrpheus_qrs` object from [mrpheus::detect_qrs()], or a
#'   named list with at least `qrs_i` (integer vector of R-peak indices) and
#'   `fs` (sampling frequency in Hz).
#'
#' @return Numeric vector of length `max(qrs$qrs_i)`. Each sample contains
#'   the instantaneous heart rate in beats per minute for the inter-beat
#'   interval it falls within. Samples before the first detected R-peak are 0.
#'
#' @seealso [mrpheus::detect_qrs()], [mrpheus::compute_hrv_sleep()]
#'
#' @export
#'
#' @examples
#' \dontrun{
#' rec <- read_philips_physlog("sub-01_physlog.log")
#' qrs <- detect_qrs(rec$C[, "v1raw"], fs = rec$HDR$sfreq)
#' hr  <- compute_hr_signal(qrs)
#'
#' plot(seq_along(hr) / rec$HDR$sfreq, hr, type = "l",
#'      xlab = "Time (s)", ylab = "HR (bpm)")
#' }
compute_hr_signal <- function(qrs) {
  if (inherits(qrs, "mrpheus_qrs")) {
    qrs_i <- qrs$qrs_i
    sfreq <- qrs$fs
  } else if (is.list(qrs) && all(c("qrs_i", "fs") %in% names(qrs))) {
    qrs_i <- qrs$qrs_i
    sfreq <- qrs$fs
  } else {
    cli::cli_abort(
      "`qrs` must be a `mrpheus_qrs` object from {.fn detect_qrs}."
    )
  }

  n <- length(qrs_i)
  if (n < 2L) {
    cli::cli_warn("Need at least two R-peaks to compute heart rate.")
    return(numeric(0))
  }

  hr <- numeric(max(qrs_i))

  for (j in seq_len(n - 1L)) {
    interval <- qrs_i[j + 1L] - qrs_i[j]
    hr[qrs_i[j]:qrs_i[j + 1L]] <- (60 * sfreq) / interval
  }

  hr
}
