#' Automatic AASM sleep staging
#'
#' Stages each 30-second epoch using a pre-trained LightGBM model originally
#' developed for YASA (Vallat & Walker, 2021) and shipped as a cross-language
#' serialised model in `inst/models/yasa_staging.txt`. Features are computed in
#' R to match the Python feature extraction pipeline exactly; bit-exact parity
#' is validated in the package test suite.
#'
#' Stages returned follow standard AASM nomenclature:
#' `W` (wake), `N1`, `N2`, `N3`, `REM`.
#'
#' @param psg An `mrpheus_psg` object from [mrpheus::prepare_psg()].
#' @param eeg_channel Character. Central EEG channel (e.g. `"EEG Fpz-Cz"`).
#'   If `NULL` (default), the first non-bad EEG channel is used.
#' @param eog_channel Character or `NULL`. EOG channel. If `NULL` (default),
#'   the first non-bad EOG channel is used (EOG features are omitted if none
#'   found).
#' @param emg_channel Character or `NULL`. Chin EMG channel. If `NULL`
#'   (default), the first non-bad EMG channel is used (EMG features are omitted
#'   if none found).
#' @param artefacts Tibble or `NULL`. Output of [mrpheus::detect_artifacts()].
#'   Artefact epochs are assigned `NA` stage and excluded from the model.
#'   If `NULL`, all epochs are staged.
#' @param model_path Character. Path to the serialised LightGBM model. Defaults
#'   to the bundled model at `system.file("models/yasa_staging.txt",
#'   package = "mrpheus")`.
#'
#' @return A tibble with one row per epoch:
#' \describe{
#'   \item{epoch}{Integer.}
#'   \item{stage}{Character. AASM stage: `W`, `N1`, `N2`, `N3`, `REM`, or
#'     `NA` (artefact).}
#'   \item{prob_W, prob_N1, prob_N2, prob_N3, prob_REM}{Numeric. Posterior
#'     class probabilities from the LightGBM model.}
#' }
#'
#' @references
#' Vallat, R., & Walker, M. P. (2021). An open-source, high-performance tool
#' for automated sleep staging. *eLife*, 10, e70092.
#' \doi{10.7554/eLife.70092}
#'
#' @seealso [mrpheus::export_hypnogram()]
#'
#' @export
stage_epochs <- function(psg,
                         eeg_channel = NULL,
                         eog_channel = NULL,
                         emg_channel = NULL,
                         artefacts   = NULL,
                         model_path  = system.file("models", "yasa_staging.txt",
                                                    package = "mrpheus")) {
  stopifnot(inherits(psg, "mrpheus_psg"))

  if (!nzchar(model_path) || !file.exists(model_path)) {
    cli::cli_abort(c(
      "LightGBM staging model not found at {.path {model_path}}.",
      "i" = "Run {.code data-raw/fetch_yasa_model.py} to download the model,",
      "i" = "then copy the output to {.path inst/models/yasa_staging.txt}."
    ))
  }

  eeg_channel <- eeg_channel %||%
    psg$channel_map$label[psg$channel_map$type == "EEG" & !psg$channel_map$bad][1]
  eog_channel <- eog_channel %||%
    psg$channel_map$label[psg$channel_map$type == "EOG" & !psg$channel_map$bad][1]
  emg_channel <- emg_channel %||%
    psg$channel_map$label[psg$channel_map$type == "EMG" & !psg$channel_map$bad][1]

  cli::cli_alert_info(
    "Staging with EEG={.val {eeg_channel}}, EOG={.val {eog_channel}}, EMG={.val {emg_channel}}"
  )

  art_epochs <- if (!is.null(artefacts)) {
    artefacts$epoch[artefacts$artefact]
  } else {
    integer(0)
  }

  features <- .extract_staging_features(psg, eeg_channel, eog_channel, emg_channel)

  model <- lightgbm::lgb.load(model_path)
  feat_mat <- as.matrix(features[, -1])

  probs <- predict(model, feat_mat, reshape = TRUE)
  colnames(probs) <- c("prob_W", "prob_N1", "prob_N2", "prob_N3", "prob_REM")
  stage_labels <- c("W", "N1", "N2", "N3", "REM")
  stage <- stage_labels[apply(probs, 1, which.max)]

  out <- tibble::tibble(
    epoch = features$epoch,
    stage = stage
  ) |>
    dplyr::bind_cols(tibble::as_tibble(probs))

  if (length(art_epochs) > 0L) {
    out$stage[out$epoch %in% art_epochs] <- NA_character_
  }

  cli::cli_alert_success("Staging complete: {nrow(out)} epochs.")
  out
}

# Feature extraction is implemented in R/staging_features.R.
