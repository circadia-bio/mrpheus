# R/preprocess_psg.R
#
# PSG preprocessing pipeline.
# Applies DC removal, notch filtering, and channel-type-specific bandpass
# filtering to the full continuous signal before re-epoching.
#
# Low-level filter functions are exported separately in signal_filters.R
# and can be used on any numeric vector independently of the pipeline.

#' Preprocess a PSG recording
#'
#' Applies the standard sleep EEG preprocessing pipeline to an `mrpheus_psg`
#' object: optional channel renaming, DC offset removal, powerline notch
#' filtering (plus harmonics), and channel-type-specific bandpass filtering.
#'
#' Filtering is applied to the full **continuous** signal stored in
#' `psg$edf$signals` before re-epoching, which avoids filter discontinuities
#' at epoch boundaries. The returned object is a new `mrpheus_psg` with the
#' same epoch length and channel map as the input.
#'
#' The individual filter steps ([remove_dc()], [detect_powerline()],
#' [notch_filter()], [bandpass_filter()]) are exported separately and can be
#' called on any numeric vector without constructing a PSG object.
#'
#' @param psg An `mrpheus_psg` object from [mrpheus::prepare_psg()].
#' @param channel_rename Named character vector mapping old channel labels to
#'   new labels, e.g. `c("C3-A2" = "C3", "O2-A1" = "O2")`. Applied before
#'   filtering; updates `channel_map` labels accordingly. Default `NULL`.
#' @param dc Logical. Remove per-channel DC offset via [remove_dc()].
#'   Default `TRUE`.
#' @param powerline_freq Integer or `NULL`. Powerline frequency in Hz (`50L`
#'   or `60L`). `NULL` (default) triggers auto-detection via
#'   [detect_powerline()] on the first clean EEG channel.
#' @param notch_harmonics Logical. Notch harmonics of `powerline_freq` up to
#'   Nyquist - 5 Hz. Default `TRUE`.
#' @param notch_bw_hz Numeric. Full bandwidth of each notch in Hz. Default `2`.
#' @param eeg_bandpass Numeric vector of length 2. Bandpass limits (Hz) applied
#'   to all EEG channels. Default `c(0.3, 35)`.
#' @param eog_bandpass Numeric vector of length 2. Bandpass limits (Hz) applied
#'   to all EOG channels. Default `c(0.3, 15)`.
#' @param emg_bandpass Numeric vector of length 2. Bandpass limits (Hz) applied
#'   to all EMG channels. Default `c(10, 99)`.
#' @param ecg_bandpass Numeric vector of length 2. Bandpass limits (Hz) applied
#'   to all ECG channels. Default `c(0.5, 40)`.
#' @param verbose Logical. Print progress messages. Default `TRUE`.
#'
#' @return A new `mrpheus_psg` object with filtered signals and re-segmented
#'   epochs. The `edf$signals` entries are replaced in-place with the filtered
#'   continuous signals so that subsequent epoch access is consistent.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' rec  <- read_edf("data/psg_001.edf")
#' psg  <- prepare_psg(rec)
#'
#' # Auto-detect powerline, apply standard filters
#' psg_clean <- preprocess_psg(psg)
#'
#' # Rename linked-ear channels to 10-20 names, fix powerline explicitly
#' psg_clean <- preprocess_psg(
#'   psg,
#'   channel_rename = c("C3-A2" = "C3", "C4-A1" = "C4",
#'                      "O1-A2" = "O1", "O2-A1" = "O2"),
#'   powerline_freq = 50L
#' )
#' }
preprocess_psg <- function(psg,
                            channel_rename  = NULL,
                            dc              = TRUE,
                            powerline_freq  = NULL,
                            notch_harmonics = TRUE,
                            notch_bw_hz     = 2,
                            eeg_bandpass    = c(0.3, 35),
                            eog_bandpass    = c(0.3, 15),
                            emg_bandpass    = c(10, 99),
                            ecg_bandpass    = c(0.5, 40),
                            verbose         = TRUE) {
  stopifnot(inherits(psg, "mrpheus_psg"))

  cmap <- psg$channel_map
  edf  <- psg$edf

  # -- Channel renaming -----------------------------------------------------------
  if (!is.null(channel_rename)) {
    old_labels <- names(channel_rename)
    found      <- old_labels %in% cmap$label

    if (any(!found)) {
      cli::cli_alert_warning(
        "channel(s) not found in channel_map, skipping: {.val {old_labels[!found]}}"
      )
    }

    for (old in old_labels[found]) {
      new <- channel_rename[[old]]
      cmap$label[cmap$label == old]                 <- new
      names(edf$signals)[names(edf$signals) == old] <- new
      edf$channels$label[edf$channels$label == old] <- new
    }

    if (verbose && any(found))
      cli::cli_alert_success("Renamed {sum(found)} channel(s)")
  }

  # -- Auto-detect powerline ------------------------------------------------------
  if (is.null(powerline_freq)) {
    eeg_ch <- cmap$label[cmap$type == "EEG" & !cmap$bad][1]

    if (is.na(eeg_ch)) {
      powerline_freq <- 50L
      if (verbose)
        cli::cli_alert_warning(
          "No clean EEG channel available for powerline detection; defaulting to 50 Hz."
        )
    } else {
      sr_ref         <- cmap$sample_rate[cmap$label == eeg_ch]
      powerline_freq <- detect_powerline(edf$signals[[eeg_ch]]$signal, sr_ref)
      if (verbose)
        cli::cli_alert_info("Detected powerline frequency: {powerline_freq} Hz")
    }
  }

  # -- Filter each non-bad channel -----------------------------------------------
  active_channels <- cmap$label[!cmap$bad]

  for (lbl in active_channels) {
    sig <- edf$signals[[lbl]]$signal
    sr  <- cmap$sample_rate[cmap$label == lbl]
    typ <- cmap$type[cmap$label == lbl]

    if (verbose)
      cli::cli_alert_info("Filtering {lbl}  [{typ} | {sr} Hz]")

    if (dc) sig <- remove_dc(sig)

    sig <- notch_filter(sig, sr,
                        freq         = powerline_freq,
                        bandwidth_hz = notch_bw_hz,
                        harmonics    = notch_harmonics)

    bp <- switch(typ,
      EEG  = eeg_bandpass,
      EOG  = eog_bandpass,
      EMG  = emg_bandpass,
      ECG  = ecg_bandpass,
      NULL
    )
    if (!is.null(bp)) sig <- bandpass_filter(sig, sr, bp[1], bp[2])

    edf$signals[[lbl]]$signal <- sig
  }

  # -- Re-epoch with filtered signals --------------------------------------------
  n_epochs <- psg$n_epochs
  epoch_s  <- psg$epoch_s

  new_epochs <- lapply(seq_len(n_epochs), function(i) {
    lapply(stats::setNames(cmap$label, cmap$label), function(lbl) {
      sig        <- edf$signals[[lbl]]$signal
      sr         <- cmap$sample_rate[cmap$label == lbl]
      ep_samples <- as.integer(epoch_s * sr)
      ep_start   <- (i - 1L) * ep_samples + 1L
      sig[ep_start:(ep_start + ep_samples - 1L)]
    })
  })

  if (verbose)
    cli::cli_alert_success(
      "Preprocessing complete -- {n_epochs} epochs ready."
    )

  structure(
    list(
      edf         = edf,
      epochs      = new_epochs,
      n_epochs    = n_epochs,
      epoch_s     = epoch_s,
      channel_map = cmap
    ),
    class = "mrpheus_psg"
  )
}
