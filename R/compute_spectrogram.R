#' Compute a time-frequency spectrogram
#'
#' Generates a short-time Fourier transform (STFT) based spectrogram for a
#' single PSG channel across all epochs. Returns power (µV²/Hz) as a matrix
#' with time (epochs) on rows and frequency on columns.
#'
#' @param psg An `mrpheus_psg` object from [prepare_psg()].
#' @param channel Character. A single channel label.
#' @param freq_range Numeric vector of length 2. Frequency range to return
#'   (Hz). Default `c(0, 40)`.
#' @param window_s Numeric. STFT window length in seconds. Default `2`.
#' @param overlap Numeric in \[0, 1). Fractional window overlap. Default `0.5`.
#' @param db Logical. Return power in decibels (`10 * log10(power)`)?
#'   Default `TRUE`.
#'
#' @return A list of class `mrpheus_spectrogram` with:
#' \describe{
#'   \item{power}{Numeric matrix. Rows = epochs, columns = frequency bins.}
#'   \item{freqs}{Numeric vector. Frequency bin centres (Hz).}
#'   \item{epochs}{Integer vector. Epoch indices.}
#'   \item{channel}{Character. The channel label.}
#'   \item{db}{Logical. Whether power is in dB.}
#' }
#'
#' @export
compute_spectrogram <- function(psg,
                                channel    = NULL,
                                freq_range = c(0, 40),
                                window_s   = 2,
                                overlap    = 0.5,
                                db         = TRUE) {
  stopifnot(inherits(psg, "mrpheus_psg"))

  if (is.null(channel)) {
    channel <- psg$channel_map$label[
      psg$channel_map$type == "EEG" & !psg$channel_map$bad
    ][1]
    cli::cli_alert_info("No channel specified; using {.val {channel}}.")
  }

  sr    <- psg$channel_map$sample_rate[psg$channel_map$label == channel]
  nfft  <- 2^ceiling(log2(window_s * sr))
  novlp <- as.integer(nfft * overlap)

  mat <- do.call(rbind, lapply(seq_along(psg$epochs), function(i) {
    sig <- psg$epochs[[i]][[channel]]
    if (is.null(sig) || length(sig) < nfft) return(rep(NA_real_, nfft / 2))
    psd <- gsignal::pwelch(sig, fs = sr, window = nfft,
                            noverlap = novlp, nfft = nfft)
    freq_idx <- psd$freq >= freq_range[1] & psd$freq <= freq_range[2]
    psd$spec[freq_idx]
  }))

  # Recompute freq vector for subsetting
  ref_psd  <- gsignal::pwelch(psg$epochs[[1]][[channel]], fs = sr,
                                window = nfft, noverlap = novlp, nfft = nfft)
  freq_idx <- ref_psd$freq >= freq_range[1] & ref_psd$freq <= freq_range[2]
  freqs    <- ref_psd$freq[freq_idx]

  if (db) mat <- 10 * log10(mat + .Machine$double.eps)

  structure(
    list(
      power   = mat,
      freqs   = freqs,
      epochs  = seq_len(nrow(mat)),
      channel = channel,
      db      = db
    ),
    class = "mrpheus_spectrogram"
  )
}

#' @export
print.mrpheus_spectrogram <- function(x, ...) {
  cli::cli_h1("mrpheus spectrogram")
  cli::cli_inform(c(
    "i" = "Channel: {.val {x$channel}}",
    "i" = "Epochs: {length(x$epochs)}",
    "i" = "Frequency range: {min(x$freqs)}-{max(x$freqs)} Hz ({length(x$freqs)} bins)",
    "i" = "Units: {if (x$db) 'dB' else 'uV2/Hz'}"
  ))
  invisible(x)
}
