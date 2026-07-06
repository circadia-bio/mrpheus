# R/compute_temporal_bandpower.R
#
# Sliding-window band power analysis across the full recording night.
# Ports calculate_temporal_bandpower() from the companion Python notebook.
#
# Each window spans `window_epochs` consecutive 30-s epochs. Windows step
# by `step_epochs` (default: window_epochs %/% 2, i.e. 50 % overlap).
# Within each window, epochs are concatenated per channel, Welch PSD is
# computed per channel, then averaged across channels before band
# integration. The result is a long-format tibble with one row per
# (window, band).

#' Compute temporal band power across a PSG recording
#'
#' Applies a sliding window over the scored epochs, concatenates the signal
#' within each window, estimates the Welch PSD averaged across channels, and
#' integrates power within standard EEG frequency bands. Returns a long-format
#' tibble suitable for plotting temporal power dynamics across the night.
#'
#' Mirrors `calculate_temporal_bandpower()` from the companion Python notebook:
#' windows step by 50 % of their length (default), Welch overlap defaults to
#' half the window length in samples (proportional to sample rate), and
#' relative power is the share of each band in the total band-summed power.
#'
#' @param psg An `mrpheus_psg` object from [mrpheus::prepare_psg()].
#' @param hypno Integer vector of length `psg$n_epochs`. Sleep stage codes
#'   following the YASA convention: `0` = Wake, `1` = N1, `2` = N2, `3` = N3,
#'   `4` = REM, `-1` = Artefact. Used to assign a dominant stage to each
#'   window.
#' @param channels Character vector. Channel labels whose PSDs are averaged.
#'   If `NULL` (default), all non-bad EEG channels are used.
#' @param window_epochs Integer. Number of 30-s epochs per window. Default
#'   `10` (5 minutes).
#' @param step_epochs Integer. Window step size in epochs. Default `NULL`,
#'   which sets the step to `window_epochs %/% 2` (50 % overlap).
#' @param bands Named list of length-2 numeric vectors (Hz). Default matches
#'   [compute_band_power()]:
#'   ```
#'   list(delta = c(0.5, 4), theta = c(4, 8), alpha = c(8, 12),
#'        sigma = c(12, 16), beta = c(16, 30), gamma = c(30, 45))
#'   ```
#' @param win_sec Numeric. Welch window length in seconds. Default `4`.
#' @param noverlap Integer or `NULL`. Welch overlap in samples. `NULL`
#'   (default) sets it to `win_sec / 2 * sample_rate`, i.e. 50 % of the
#'   Welch window — matching the notebook's `noverlap = int(2 * sfreq)`.
#' @param nfft Integer or `NULL`. FFT length. `NULL` (default) uses the next
#'   power of 2 above `win_sec * sample_rate`.
#'
#' @return A tibble with columns:
#' \describe{
#'   \item{time_hours}{Numeric. Start time of the window in hours from
#'     recording onset.}
#'   \item{epoch_start}{Integer. First epoch index in the window (1-based).}
#'   \item{epoch_end}{Integer. Last epoch index in the window (1-based).}
#'   \item{dominant_stage}{Integer. Most common stage code across the window
#'     (artefact epochs excluded). `NA` if all epochs are artefacts.}
#'   \item{band}{Character. Band name.}
#'   \item{power}{Numeric. Mean band power across channels (V^2/Hz).}
#'   \item{relative_power}{Numeric. Band power as a fraction of total
#'     band-summed power at this time point.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' rec   <- read_edf("SC4001E0-PSG.edf")
#' psg   <- prepare_psg(rec) |> preprocess_psg()
#' hypno <- stage_epochs(psg)$stage
#'
#' tbp <- compute_temporal_bandpower(psg, hypno)
#'
#' # Plot delta power across the night
#' library(ggplot2)
#' tbp |>
#'   dplyr::filter(band == "delta") |>
#'   ggplot(aes(time_hours, relative_power)) +
#'   geom_line() +
#'   labs(x = "Time (hours)", y = "Relative delta power")
#' }
compute_temporal_bandpower <- function(psg,
                                        hypno,
                                        channels      = NULL,
                                        window_epochs = 10L,
                                        step_epochs   = NULL,
                                        bands         = list(
                                          delta = c(0.5,  4),
                                          theta = c(4,    8),
                                          alpha = c(8,   12),
                                          sigma = c(12,  16),
                                          beta  = c(16,  30),
                                          gamma = c(30,  45)
                                        ),
                                        win_sec  = 4,
                                        noverlap = NULL,
                                        nfft     = NULL) {
  stopifnot(inherits(psg, "mrpheus_psg"))

  if (length(hypno) != psg$n_epochs) {
    cli::cli_abort(
      "Length of {.arg hypno} ({length(hypno)}) must equal \\
       {.code psg$n_epochs} ({psg$n_epochs})."
    )
  }

  if (is.null(channels)) {
    channels <- psg$channel_map$label[
      psg$channel_map$type == "EEG" & !psg$channel_map$bad
    ]
  }
  if (length(channels) == 0L)
    cli::cli_abort("No valid EEG channels found.")

  if (is.null(step_epochs))
    step_epochs <- max(1L, as.integer(window_epochs) %/% 2L)

  n_epochs <- psg$n_epochs
  epoch_s  <- psg$epoch_s

  if (as.integer(window_epochs) > n_epochs)
    cli::cli_abort(
      "{.arg window_epochs} ({window_epochs}) exceeds total epochs ({n_epochs})."
    )

  starts <- seq(1L, n_epochs - as.integer(window_epochs) + 1L,
                by = as.integer(step_epochs))

  if (length(starts) == 0L)
    cli::cli_abort(
      "{.arg window_epochs} ({window_epochs}) exceeds total epochs ({n_epochs})."
    )

  # Per-channel sample rate (used to compute noverlap / nfft defaults)
  sr_map <- stats::setNames(psg$channel_map$sample_rate, psg$channel_map$label)

  rows <- lapply(starts, function(s) {
    e_end <- min(s + as.integer(window_epochs) - 1L, n_epochs)

    # Time in hours at window start (0-indexed epoch count, matching notebook)
    time_h <- ((s - 1L) * epoch_s) / 3600

    # Dominant stage: most common among non-artefact epochs in window
    win_hypno    <- hypno[s:e_end]
    valid_stages <- win_hypno[win_hypno >= 0L]
    dom_stage    <- if (length(valid_stages) > 0L) {
      as.integer(names(sort(table(valid_stages), decreasing = TRUE))[1L])
    } else {
      NA_integer_
    }

    # Average PSD across channels
    psd_list <- lapply(channels, function(ch) {
      sr      <- sr_map[[ch]]
      nov     <- if (!is.null(noverlap)) noverlap else as.integer(win_sec / 2 * sr)
      nf      <- if (!is.null(nfft))    nfft      else
                   2L ^ ceiling(log2(as.integer(win_sec * sr)))

      sig_win <- unlist(lapply(s:e_end, function(i) psg$epochs[[i]][[ch]]),
                        use.names = FALSE)
      .welch_psd_bp(sig_win, fs = sr,
                    win_sec  = win_sec,
                    noverlap = nov,
                    nfft     = nf)
    })

    avg_spec <- Reduce(`+`, lapply(psd_list, `[[`, "spec")) / length(psd_list)
    freq     <- psd_list[[1L]]$freq

    # Integrate per band
    band_power <- vapply(names(bands), function(bn) {
      b   <- bands[[bn]]
      idx <- freq >= b[1L] & freq < b[2L]
      if (!any(idx)) return(NA_real_)
      pracma::trapz(freq[idx], avg_spec[idx])
    }, numeric(1L))

    total <- sum(band_power, na.rm = TRUE)

    tibble::tibble(
      time_hours     = time_h,
      epoch_start    = s,
      epoch_end      = e_end,
      dominant_stage = dom_stage,
      band           = names(bands),
      power          = band_power,
      relative_power = band_power / total
    )
  })

  dplyr::bind_rows(rows)
}
