#' Compute EEG band power per epoch
#'
#' Estimates power spectral density (PSD) using Welch's method and integrates
#' power within standard EEG frequency bands (delta, theta, alpha, sigma, beta,
#' gamma) for each epoch and each specified channel. Mirrors the band-power
#' feature extraction used by YASA's staging pipeline.
#'
#' @param psg An `mrpheus_psg` object from [mrpheus::prepare_psg()].
#' @param channels Character vector. EEG channel labels. If `NULL` (default),
#'   all non-bad EEG channels are used.
#' @param bands Named list of length-2 numeric vectors defining frequency bands.
#'   Default:
#'   ```
#'   list(delta = c(0.5, 4), theta = c(4, 8), alpha = c(8, 13),
#'        sigma = c(13, 16), beta = c(16, 30), gamma = c(30, 40))
#'   ```
#' @param relative Logical. If `TRUE`, return relative band power (band / total
#'   power in 0.5–40 Hz). Default `FALSE`.
#' @param window_s Numeric. Welch window length in seconds. Default `4`.
#' @param overlap Numeric in \[0, 1). Fractional overlap between Welch windows.
#'   Default `0.5`.
#'
#' @return A tibble with columns `epoch`, `channel`, one column per band, and
#'   `total_power`. Units are µV²/Hz (or dimensionless if `relative = TRUE`).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' bp <- compute_band_power(psg)
#' bp <- compute_band_power(psg, channels = "EEG Fpz-Cz", relative = TRUE)
#' }
compute_band_power <- function(psg,
                               channels  = NULL,
                               bands     = list(
                                 delta = c(0.5, 4),
                                 theta = c(4, 8),
                                 alpha = c(8, 13),
                                 sigma = c(13, 16),
                                 beta  = c(16, 30),
                                 gamma = c(30, 40)
                               ),
                               relative  = FALSE,
                               window_s  = 4,
                               overlap   = 0.5) {
  stopifnot(inherits(psg, "mrpheus_psg"))

  if (is.null(channels)) {
    channels <- psg$channel_map$label[
      psg$channel_map$type == "EEG" & !psg$channel_map$bad
    ]
  }

  rows <- lapply(seq_along(psg$epochs), function(i) {
    ep <- psg$epochs[[i]]

    lapply(channels, function(ch) {
      sig <- ep[[ch]]
      sr  <- psg$channel_map$sample_rate[psg$channel_map$label == ch]
      if (is.null(sig) || length(sig) < 2) return(NULL)

      nfft   <- 2^ceiling(log2(window_s * sr))
      novlap <- as.integer(nfft * overlap)

      psd <- gsignal::pwelch(sig, fs = sr, window = nfft,
                              noverlap = novlap, nfft = nfft)

      band_power <- vapply(bands, function(b) {
        idx <- psd$freq >= b[1] & psd$freq < b[2]
        if (!any(idx)) return(NA_real_)
        pracma::trapz(psd$freq[idx], psd$spec[idx])
      }, numeric(1))

      total_range <- psd$freq >= 0.5 & psd$freq <= 40
      total       <- pracma::trapz(psd$freq[total_range],
                                    psd$spec[total_range])

      if (relative) band_power <- band_power / total

      row <- as.list(band_power)
      row$epoch       <- i
      row$channel     <- ch
      row$total_power <- total
      tibble::as_tibble(row)
    })
  })

  dplyr::bind_rows(unlist(rows, recursive = FALSE)) |>
    dplyr::relocate(epoch, channel)
}
