#' Detect artefact epochs in a PSG recording
#'
#' Flags epochs containing likely artefacts based on amplitude thresholds,
#' excessive high-frequency (muscle) power, and movement contamination. Operates
#' on the EEG channels by default. Returns a logical vector (one value per epoch)
#' and optionally attaches per-epoch artefact metrics for inspection.
#'
#' @param psg An `mrpheus_psg` object from [prepare_psg()].
#' @param channels Character vector. Channel labels to evaluate. If `NULL`
#'   (default), all non-bad EEG channels are used.
#' @param amp_threshold_uv Numeric. Peak-to-peak amplitude threshold in µV.
#'   Epochs exceeding this in any selected channel are flagged. Default `500`.
#' @param hf_band Numeric vector of length 2. Frequency band (Hz) considered
#'   high-frequency / muscle contamination. Default `c(40, 100)`.
#' @param hf_percentile Numeric. Epochs whose HF power exceeds this percentile
#'   (across all epochs) are flagged. Default `0.99`.
#' @param verbose Logical. Print summary. Default `TRUE`.
#'
#' @return A tibble with one row per epoch and columns:
#' \describe{
#'   \item{epoch}{Integer. Epoch index (1-based).}
#'   \item{artefact}{Logical. `TRUE` if the epoch is flagged as artefact.}
#'   \item{reason}{Character. Comma-separated reasons for flagging, or `NA`.}
#'   \item{peak_to_peak_uv}{Numeric. Maximum peak-to-peak amplitude across
#'     selected channels.}
#'   \item{hf_power_db}{Numeric. Mean HF band power (dB) across selected
#'     channels.}
#' }
#'
#' @export
detect_artifacts <- function(psg,
                             channels          = NULL,
                             amp_threshold_uv  = 500,
                             hf_band           = c(40, 100),
                             hf_percentile     = 0.99,
                             verbose           = TRUE) {
  stopifnot(inherits(psg, "mrpheus_psg"))

  # Resolve channels
  if (is.null(channels)) {
    channels <- psg$channel_map$label[
      psg$channel_map$type == "EEG" & !psg$channel_map$bad
    ]
  }
  if (length(channels) == 0L) {
    cli::cli_abort("No valid EEG channels available for artefact detection.")
  }

  sr <- psg$channel_map$sample_rate[psg$channel_map$label == channels[1]]

  results <- lapply(seq_along(psg$epochs), function(i) {
    ep <- psg$epochs[[i]]

    # Amplitude check
    pp <- max(vapply(channels, function(ch) {
      s <- ep[[ch]]
      if (is.null(s) || all(is.na(s))) return(0)
      diff(range(s, na.rm = TRUE))
    }, numeric(1)))

    # HF power check
    hf_pwr <- mean(vapply(channels, function(ch) {
      s <- ep[[ch]]
      if (is.null(s) || length(s) < 2) return(0)
      psd <- gsignal::pwelch(s, fs = sr, nfft = min(length(s), 256L),
                             window = min(length(s), 256L))
      freq_idx <- psd$freq >= hf_band[1] & psd$freq <= hf_band[2]
      if (!any(freq_idx)) return(0)
      mean(10 * log10(psd$spec[freq_idx] + .Machine$double.eps))
    }, numeric(1)))

    list(epoch = i, pp = pp, hf = hf_pwr)
  })

  pp_vals  <- vapply(results, `[[`, numeric(1), "pp")
  hf_vals  <- vapply(results, `[[`, numeric(1), "hf")
  hf_cutoff <- quantile(hf_vals, hf_percentile, na.rm = TRUE)

  artefact <- pp_vals > amp_threshold_uv | hf_vals > hf_cutoff

  reason <- mapply(function(pp, hf) {
    r <- c(
      if (pp > amp_threshold_uv) "amplitude" else NULL,
      if (hf > hf_cutoff)        "high_freq"  else NULL
    )
    if (length(r) == 0L) NA_character_ else paste(r, collapse = ",")
  }, pp_vals, hf_vals, SIMPLIFY = TRUE)

  out <- tibble::tibble(
    epoch           = seq_along(artefact),
    artefact        = artefact,
    reason          = reason,
    peak_to_peak_uv = pp_vals,
    hf_power_db     = hf_vals
  )

  if (verbose) {
    n_art <- sum(artefact)
    cli::cli_alert_info(
      "Artefact epochs: {n_art} / {nrow(out)} ({round(100 * n_art / nrow(out), 1)}%)"
    )
  }

  out
}
