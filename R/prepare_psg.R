#' Prepare a PSG recording for analysis
#'
#' Takes an `mrpheus_edf` object and segments it into standard epochs,
#' performs a channel inventory (classifying signals by type), and flags
#' channels that appear flat or likely bad. This is the standard entry point
#' before any downstream analyses (spectral, event detection, staging).
#'
#' @param edf An `mrpheus_edf` object from [mrpheus::read_edf()].
#' @param epoch_s Numeric. Epoch length in seconds. Default `30` (standard
#'   AASM epoch).
#' @param eeg_pattern Character. Regex pattern to identify EEG channels.
#'   Default `"EEG|C3|C4|F3|F4|O1|O2|Fpz|Pz"`.
#' @param eog_pattern Character. Regex pattern to identify EOG channels.
#'   Default `"EOG|ROC|LOC"`.
#' @param emg_pattern Character. Regex pattern to identify EMG channels.
#'   Default `"EMG|chin|Chin"`.
#' @param ecg_pattern Character. Regex pattern to identify ECG/EKG channels.
#'   Default `"ECG|EKG"`.
#' @param resp_pattern Character. Regex pattern to identify respiratory channels.
#'   Default `"Thor|Abdo|Flow|SpO2|airflow"`.
#' @param flat_threshold Numeric. Variance below this value flags a channel as
#'   flat/bad. Default `1e-6`.
#'
#' @return A list of class `mrpheus_psg` with components:
#' \describe{
#'   \item{edf}{The original `mrpheus_edf` object.}
#'   \item{epochs}{List. Each element is one epoch (30 s by default), itself a
#'     named list of channel vectors.}
#'   \item{n_epochs}{Integer. Total number of complete epochs.}
#'   \item{epoch_s}{Numeric. Epoch duration in seconds.}
#'   \item{channel_map}{Data frame. Channel label, detected type
#'     (EEG/EOG/EMG/ECG/RESP/OTHER), sample rate, and `bad` flag.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' rec  <- read_edf("data/psg_001.edf")
#' psg  <- prepare_psg(rec)
#' psg$channel_map
#' }
prepare_psg <- function(edf,
                        epoch_s        = 30,
                        eeg_pattern    = "EEG|C3|C4|F3|F4|O1|O2|Fpz|Pz",
                        eog_pattern    = "EOG|ROC|LOC",
                        emg_pattern    = "EMG|chin|Chin",
                        ecg_pattern    = "ECG|EKG",
                        resp_pattern   = "Thor|Abdo|Flow|SpO2|airflow",
                        flat_threshold = 1e-6) {
  stopifnot(inherits(edf, "mrpheus_edf"))

  labels <- edf$channels$label
  type <- dplyr::case_when(
    grepl(eeg_pattern,  labels, ignore.case = FALSE) ~ "EEG",
    grepl(eog_pattern,  labels, ignore.case = FALSE) ~ "EOG",
    grepl(emg_pattern,  labels, ignore.case = FALSE) ~ "EMG",
    grepl(ecg_pattern,  labels, ignore.case = FALSE) ~ "ECG",
    grepl(resp_pattern, labels, ignore.case = FALSE) ~ "RESP",
    TRUE                                              ~ "OTHER"
  )

  bad <- vapply(labels, function(lbl) {
    sig <- edf$signals[[lbl]]$signal
    if (is.null(sig) || length(sig) < 2) return(TRUE)
    var(sig, na.rm = TRUE) < flat_threshold
  }, logical(1))

  if (any(bad)) {
    cli::cli_alert_warning(
      "Flat/bad channels detected: {.val {labels[bad]}}"
    )
  }

  channel_map <- tibble::tibble(
    label       = labels,
    type        = type,
    sample_rate = edf$channels$sample_rate,
    bad         = bad
  )

  ref_label <- labels[type == "EEG" & !bad][1]
  if (is.na(ref_label)) {
    cli::cli_abort("No valid EEG channel found for epoch segmentation.")
  }

  ref_sig           <- edf$signals[[ref_label]]$signal
  ref_sr            <- channel_map$sample_rate[channel_map$label == ref_label]
  samples_per_epoch <- epoch_s * ref_sr
  n_epochs          <- floor(length(ref_sig) / samples_per_epoch)

  cli::cli_alert_info(
    "Segmenting into {n_epochs} epochs of {epoch_s} s each."
  )

  epochs <- lapply(seq_len(n_epochs), function(i) {
    lapply(stats::setNames(labels, labels), function(lbl) {
      sig        <- edf$signals[[lbl]]$signal
      sr         <- channel_map$sample_rate[channel_map$label == lbl]
      ep_samples <- epoch_s * sr
      ep_start   <- (i - 1L) * ep_samples + 1L
      sig[ep_start:(ep_start + ep_samples - 1L)]
    })
  })

  structure(
    list(
      edf         = edf,
      epochs      = epochs,
      n_epochs    = n_epochs,
      epoch_s     = epoch_s,
      channel_map = channel_map
    ),
    class = "mrpheus_psg"
  )
}

#' @export
print.mrpheus_psg <- function(x, ...) {
  cli::cli_h1("mrpheus PSG")
  cli::cli_inform(c(
    "i" = "Epochs: {x$n_epochs} x {x$epoch_s} s",
    "i" = "Recording: {round(x$n_epochs * x$epoch_s / 3600, 2)} h"
  ))
  print(x$channel_map, ...)
  invisible(x)
}
