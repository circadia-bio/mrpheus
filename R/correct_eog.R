# R/correct_eog.R
#
# EOG artefact correction for PSG recordings.
#
# Two complementary methods are provided:
#
#   correct_eog_regression() -- multiple linear regression of each EEG
#     channel on the EOG reference channels. Equivalent to MNE's SSP
#     projection approach but requires no blink event detection: the full
#     continuous EOG signal acts as the regressor, which removes both
#     transient blinks and slow ocular drift.
#
#   correct_eog_ica() -- Independent Component Analysis via fastICA.
#     Components whose time courses correlate with the EOG above a
#     threshold are identified and their contributions subtracted from
#     the EEG. Mirrors MNE's ICA pipeline (threshold = 0.35,
#     n_components = min(6, n_eeg) by default).
#
# Both functions operate on the full continuous signal stored in
# psg$edf$signals and return a new mrpheus_psg with re-segmented epochs.

# ── Internal re-epoch helper ──────────────────────────────────────────────────

.rebuild_epochs <- function(edf, cmap, n_epochs, epoch_s) {
  lapply(seq_len(n_epochs), function(i) {
    lapply(stats::setNames(cmap$label, cmap$label), function(lbl) {
      sig        <- edf$signals[[lbl]]$signal
      sr         <- cmap$sample_rate[cmap$label == lbl]
      ep_samples <- as.integer(epoch_s * sr)
      ep_start   <- (i - 1L) * ep_samples + 1L
      sig[ep_start:(ep_start + ep_samples - 1L)]
    })
  })
}

# ── correct_eog_regression ────────────────────────────────────────────────────

#' Remove EOG artefacts using linear regression
#'
#' For each EEG channel, fits a multiple linear regression on all EOG reference
#' channels and subtracts the fitted EOG component, leaving the residual as the
#' cleaned signal. This removes both transient blink artefacts and slow ocular
#' drift without requiring blink event detection.
#'
#' This is the pure-R equivalent of MNE's `compute_proj_eog()` / SSP
#' projection pipeline. Both approaches project out the subspace spanned by
#' the EOG channels; regression does so directly on the continuous signal.
#'
#' @param psg An `mrpheus_psg` object from [mrpheus::prepare_psg()].
#' @param eog_channels Character vector. EOG channel labels to use as
#'   regressors. If `NULL` (default), all channels of type `"EOG"` in
#'   `psg$channel_map` are used.
#' @param eeg_channels Character vector. EEG channel labels to clean. If
#'   `NULL` (default), all non-bad `"EEG"` channels are used.
#' @param verbose Logical. Print progress messages. Default `TRUE`.
#'
#' @return A new `mrpheus_psg` with cleaned EEG signals and re-segmented
#'   epochs. EOG and other channels are unchanged.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' rec   <- read_edf("psg.edf")
#' psg   <- prepare_psg(rec) |> preprocess_psg()
#' clean <- correct_eog_regression(psg)
#' }
correct_eog_regression <- function(psg,
                                    eog_channels = NULL,
                                    eeg_channels = NULL,
                                    verbose      = TRUE) {
  stopifnot(inherits(psg, "mrpheus_psg"))

  cmap <- psg$channel_map
  edf  <- psg$edf

  if (is.null(eog_channels))
    eog_channels <- cmap$label[cmap$type == "EOG" & !cmap$bad]
  if (is.null(eeg_channels))
    eeg_channels <- cmap$label[cmap$type == "EEG" & !cmap$bad]

  if (length(eog_channels) == 0L)
    cli::cli_abort("No EOG channels found. Supply {.arg eog_channels} explicitly.")
  if (length(eeg_channels) == 0L)
    cli::cli_abort("No EEG channels found. Supply {.arg eeg_channels} explicitly.")

  if (verbose)
    cli::cli_alert_info(
      "Regression EOG correction: {length(eog_channels)} EOG reference(s), \\
       {length(eeg_channels)} EEG channel(s)."
    )

  # Design matrix: intercept + EOG channels (continuous signal)
  eog_mat <- do.call(cbind, lapply(eog_channels, function(ch)
    edf$signals[[ch]]$signal))
  X <- cbind(1, eog_mat)

  for (ch in eeg_channels) {
    eeg_sig                  <- edf$signals[[ch]]$signal
    beta                     <- stats::.lm.fit(X, eeg_sig)$coefficients
    edf$signals[[ch]]$signal <- eeg_sig - X %*% beta
  }

  if (verbose)
    cli::cli_alert_success("Regression correction applied.")

  structure(
    list(
      edf         = edf,
      epochs      = .rebuild_epochs(edf, cmap, psg$n_epochs, psg$epoch_s),
      n_epochs    = psg$n_epochs,
      epoch_s     = psg$epoch_s,
      channel_map = cmap
    ),
    class = "mrpheus_psg"
  )
}

# ── correct_eog_ica ───────────────────────────────────────────────────────────

#' Remove EOG artefacts using ICA
#'
#' Decomposes the EEG into independent components via FastICA, identifies
#' components whose time courses correlate with the EOG channels above
#' `threshold`, subtracts their contribution from the EEG, and returns a
#' cleaned `mrpheus_psg`. Mirrors MNE's ICA EOG rejection pipeline with
#' default settings matching the companion Python notebook
#' (`n_components = min(6, n_eeg)`, `threshold = 0.35`).
#'
#' **Note:** FastICA (this implementation) and MNE's default InfoMax ICA find
#' different decompositions; the identified components and final signal will
#' therefore differ from MNE output, but artifact removal performance is
#' comparable.
#'
#' @param psg An `mrpheus_psg` object from [mrpheus::prepare_psg()].
#' @param eog_channels Character vector. EOG channel labels used as the
#'   artifact reference. If `NULL` (default), all `"EOG"` channels are used.
#' @param eeg_channels Character vector. EEG channels to decompose and clean.
#'   If `NULL` (default), all non-bad `"EEG"` channels are used.
#' @param n_components Integer or `NULL`. Number of ICA components. `NULL`
#'   (default) uses `min(6L, n_eeg_channels)`, matching the notebook.
#' @param threshold Numeric. Absolute Pearson correlation threshold above
#'   which a component is flagged as EOG-related. Default `0.35`.
#' @param fun Character. Contrast function passed to [fastICA::fastICA()].
#'   `"logcosh"` (default) or `"exp"`. The R method is used internally for
#'   robustness across channel counts; this is slower than the C backend but
#'   avoids matrix-conformality errors on small channel sets.
#' @param verbose Logical. Print progress messages. Default `TRUE`.
#'
#' @return A new `mrpheus_psg` with cleaned EEG signals and re-segmented
#'   epochs. EOG and other channels are unchanged.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' rec   <- read_edf("psg.edf")
#' psg   <- prepare_psg(rec) |> preprocess_psg()
#' clean <- correct_eog_ica(psg)
#'
#' # Stricter threshold
#' clean <- correct_eog_ica(psg, threshold = 0.5)
#' }
correct_eog_ica <- function(psg,
                              eog_channels = NULL,
                              eeg_channels = NULL,
                              n_components = NULL,
                              threshold    = 0.35,
                              fun          = "logcosh",
                              verbose      = TRUE) {
  stopifnot(inherits(psg, "mrpheus_psg"))

  cmap <- psg$channel_map
  edf  <- psg$edf

  if (is.null(eog_channels))
    eog_channels <- cmap$label[cmap$type == "EOG" & !cmap$bad]
  if (is.null(eeg_channels))
    eeg_channels <- cmap$label[cmap$type == "EEG" & !cmap$bad]

  if (length(eog_channels) == 0L)
    cli::cli_abort("No EOG channels found. Supply {.arg eog_channels} explicitly.")
  if (length(eeg_channels) == 0L)
    cli::cli_abort("No EEG channels found. Supply {.arg eeg_channels} explicitly.")

  n_eeg        <- length(eeg_channels)
  if (n_eeg < 2L)
    cli::cli_abort(
      "ICA requires at least 2 EEG channels; found {n_eeg}. \n  Use {.fn correct_eog_regression} for single-channel recordings."
    )

  n_components <- if (is.null(n_components)) min(6L, n_eeg) else
                    min(as.integer(n_components), n_eeg - 1L)

  if (verbose)
    cli::cli_alert_info(
      "ICA EOG correction: {n_components} component(s), \\
       threshold = {threshold}."
    )

  # EEG data matrix: samples x channels (always a matrix, even for 1 channel)
  eeg_mat <- do.call(cbind, lapply(eeg_channels, function(ch)
    edf$signals[[ch]]$signal))
  if (!is.matrix(eeg_mat)) eeg_mat <- matrix(eeg_mat, ncol = 1L)

  ica <- fastICA::fastICA(eeg_mat,
                           n.comp   = n_components,
                           alg.typ  = "parallel",
                           fun      = fun,
                           method   = "R",
                           verbose  = FALSE)

  # EOG reference matrix: samples x eog_channels
  eog_mat <- do.call(cbind, lapply(eog_channels, function(ch)
    edf$signals[[ch]]$signal))

  # Max absolute correlation of each IC with any EOG channel
  cors <- apply(ica$S, 2L, function(comp) {
    if (stats::sd(comp) < .Machine$double.eps) return(0)
    max(abs(stats::cor(comp, eog_mat)))
  })

  bad <- which(cors > threshold)

  if (length(bad) == 0L) {
    if (verbose)
      cli::cli_alert_warning(
        "No ICA components exceeded the correlation threshold ({threshold}). \\
         Returning data unchanged."
      )
    return(structure(
      list(edf = edf, epochs = psg$epochs,
           n_epochs = psg$n_epochs, epoch_s = psg$epoch_s, channel_map = cmap),
      class = "mrpheus_psg"
    ))
  }

  if (verbose)
    cli::cli_alert_info(
      "Removing {length(bad)} component(s) \\
       (max EOG r = {round(max(cors[bad]), 3)})."
    )

  # ICA model: eeg_mat = ica$S %*% t(ica$A)
  # Remove bad components: subtract their contribution
  S_bad      <- ica$S[, bad, drop = FALSE]   # samples x n_bad
  A_bad      <- ica$A[, bad, drop = FALSE]   # n_eeg   x n_bad
  correction <- S_bad %*% t(A_bad)           # samples x n_eeg

  for (i in seq_along(eeg_channels)) {
    edf$signals[[eeg_channels[i]]]$signal <-
      edf$signals[[eeg_channels[i]]]$signal - correction[, i]
  }

  if (verbose)
    cli::cli_alert_success("ICA correction applied.")

  structure(
    list(
      edf         = edf,
      epochs      = .rebuild_epochs(edf, cmap, psg$n_epochs, psg$epoch_s),
      n_epochs    = psg$n_epochs,
      epoch_s     = psg$epoch_s,
      channel_map = cmap
    ),
    class = "mrpheus_psg"
  )
}
