#' Automatic AASM sleep staging
#'
#' Stages each 30-second epoch with the pre-trained LightGBM classifiers that
#' ship with YASA (Vallat & Walker, 2021), converted to LightGBM text files in
#' `inst/models/`. As in YASA, the classifier is chosen from the channels that
#' are supplied: EEG only, EEG + EOG, EEG + EMG or EEG + EOG + EMG, each with or
#' without the demographic features (`age`, `male`). Features are computed in R
#' to match the Python feature extraction and are rounded to single precision
#' before prediction, as YASA does. See the `NEWS.md` entry for the agreement
#' that was measured against YASA and for its limits.
#'
#' Stages returned follow standard AASM nomenclature:
#' `W` (wake), `N1`, `N2`, `N3`, `REM`.
#'
#' @param psg An `mrpheus_psg` object from [mrpheus::prepare_psg()].
#' @param eeg_channel Character. Central EEG channel (e.g. `"EEG Fpz-Cz"`).
#'   If `NULL` (default), the first non-bad EEG channel is used.
#' @param eog_channel Character, `NA` or `NULL`. EOG channel. If `NULL`
#'   (default), the first non-bad EOG channel is used. Without an EOG channel
#'   (`NA`, or none found) the EEG-only or EEG + EMG classifier is used.
#' @param emg_channel Character, `NA` or `NULL`. Chin EMG channel. If `NULL`
#'   (default), the first non-bad EMG channel is used. Without an EMG channel
#'   (`NA`, or none found) the classifier without EMG is used.
#' @param artefacts Tibble or `NULL`. Output of [mrpheus::detect_artifacts()].
#'   Artefact epochs are assigned `NA` stage. The classifier still sees every
#'   epoch, so the rolling-window features of neighbouring epochs are
#'   unchanged. If `NULL`, all epochs are staged.
#' @param model_path Character or `NULL`. Path to a serialised LightGBM model
#'   (plain text or `.gz`), for example one of your own. `NULL` (default) picks
#'   the bundled YASA classifier that matches the channels and `age` / `male`.
#' @param age Numeric or `NULL`. Age in years, between 0 and 120. Together with
#'   `male` it selects the `+demo` classifiers. Truncated to an integer, as in
#'   YASA.
#' @param male Logical, `0` or `1`, or `NULL`. Sex (`1` / `TRUE` for male).
#'   Must be supplied together with `age`.
#'
#' @return A tibble with one row per epoch:
#' \describe{
#'   \item{epoch}{Integer.}
#'   \item{stage}{Character. AASM stage: `W`, `N1`, `N2`, `N3`, `REM`, or
#'     `NA` (artefact).}
#'   \item{prob_W, prob_N1, prob_N2, prob_N3, prob_REM}{Numeric. Posterior
#'     class probabilities from the LightGBM model.}
#' }
#' The classifier that was used is attached as the attribute `"classifier"`.
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
                         model_path  = NULL,
                         age         = NULL,
                         male        = NULL) {
  stopifnot(inherits(psg, "mrpheus_psg"))

  # -- Demographic metadata (YASA's "+demo" classifiers need both) ------------
  if (xor(is.null(age), is.null(male))) {
    cli::cli_abort("Supply both {.arg age} and {.arg male}, or neither.")
  }
  has_demo <- !is.null(age)
  if (has_demo) {
    if (!is.numeric(age) || length(age) != 1L || !(age > 0 && age < 120))
      cli::cli_abort("{.arg age} must be a single number between 0 and 120.")
    male <- as.integer(male)
    if (length(male) != 1L || !isTRUE(male %in% 0:1))
      cli::cli_abort("{.arg male} must be 0 or 1 (or FALSE / TRUE).")
    age <- as.integer(age)               # YASA: features[c].astype(int)
  }

  # -- Channels -----------------------------------------------------------------
  first_good <- function(type) {
    psg$channel_map$label[psg$channel_map$type == type & !psg$channel_map$bad][1]
  }
  eeg_channel <- eeg_channel %||% first_good("EEG")
  eog_channel <- eog_channel %||% first_good("EOG")
  emg_channel <- emg_channel %||% first_good("EMG")
  if (is.na(eeg_channel))
    cli::cli_abort("No EEG channel available. Check `prepare_psg()` channel map.")

  has_eog <- !is.na(eog_channel)
  has_emg <- !is.na(emg_channel)

  # -- Classifier ---------------------------------------------------------------
  key <- .yasa_model_key(has_eog, has_emg, has_demo)
  classifier <- paste0("clf_", gsub("_", "+", key))
  if (is.null(model_path)) {
    model_path <- system.file("models", paste0("yasa_", key, ".txt.gz"),
                              package = "mrpheus")
  } else {
    classifier <- basename(model_path)
  }
  if (!nzchar(model_path) || !file.exists(model_path)) {
    cli::cli_abort(c(
      "LightGBM staging model not found: {.path {model_path}}.",
      "i" = "Run {.code data-raw/fetch_yasa_models.py} and copy its output into {.path inst/models/}."
    ))
  }

  cli::cli_alert_info(paste0(
    "Staging with EEG={.val {eeg_channel}}, EOG={.val {eog_channel}}, ",
    "EMG={.val {emg_channel}} (YASA classifier {.val {classifier}})"
  ))

  art_epochs <- if (!is.null(artefacts)) {
    artefacts$epoch[artefacts$artefact]
  } else {
    integer(0)
  }

  # -- Features -----------------------------------------------------------------
  features <- .extract_staging_features(psg, eeg_channel, eog_channel, emg_channel)
  feat_mat <- as.matrix(features[, -1])
  if (has_demo) {
    feat_mat <- cbind(feat_mat, age = age, male = male)
  }

  # LightGBM R passes matrices by position, not by name: use the column order
  # stored in the model file.
  model_lines <- .read_model_lines(model_path)
  feat_line   <- grep("^feature_names=", model_lines, value = TRUE)
  if (length(feat_line) != 1L)
    cli::cli_abort("Could not read {.field feature_names} from {.path {model_path}}.")
  model_feat_order <- strsplit(sub("^feature_names=", "", feat_line), " ")[[1L]]
  absent <- setdiff(model_feat_order, colnames(feat_mat))
  if (length(absent) > 0L) {
    cli::cli_abort(c(
      "The classifier needs feature(s) that were not computed: {.val {absent}}.",
      "i" = "For the {.val +demo} classifiers supply both {.arg age} and {.arg male}."
    ))
  }
  feat_mat <- feat_mat[, model_feat_order, drop = FALSE]

  # YASA rounds all float64 features to float32 before predicting. Integer
  # columns (zero-crossing counts, age, male) are exactly representable.
  feat_mat[] <- .as_float32(feat_mat)

  model <- lightgbm::lgb.load(model_str = paste(model_lines, collapse = "\n"))
  probs <- predict(model, feat_mat)
  # LightGBM v4.x returns a matrix directly for multiclass (no reshape needed).
  # Older versions returned a flat vector that needed matrix(..., byrow=TRUE).
  if (!is.matrix(probs))
    probs <- matrix(probs, ncol = 5L, byrow = TRUE)
  # The classifiers were trained with sklearn's LabelEncoder, which orders the
  # classes alphabetically: N1=0, N2=1, N3=2, R=3, W=4.
  colnames(probs) <- c("prob_N1", "prob_N2", "prob_N3", "prob_REM", "prob_W")
  stage_labels <- c("N1", "N2", "N3", "REM", "W")
  stage <- stage_labels[apply(probs, 1, which.max)]

  out <- tibble::tibble(
    epoch = features$epoch,
    stage = stage
  ) |>
    dplyr::bind_cols(tibble::as_tibble(probs))

  if (length(art_epochs) > 0L) {
    out$stage[out$epoch %in% art_epochs] <- NA_character_
  }

  attr(out, "classifier") <- classifier
  cli::cli_alert_success("Staging complete: {nrow(out)} epochs.")
  out
}

# ── Internal helpers ─────────────────────────────────────────────────────────

# Bundled classifier key, e.g. "eeg_eog_emg_demo" (YASA: clf_eeg+eog+emg+demo).
.yasa_model_key <- function(has_eog, has_emg, has_demo) {
  paste0("eeg",
         if (has_eog)  "_eog",
         if (has_emg)  "_emg",
         if (has_demo) "_demo")
}

# Lines of a LightGBM text model, plain or gzip-compressed.
.read_model_lines <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  readLines(con, warn = FALSE)
}

# Round doubles to the nearest IEEE single and return them as doubles, as
# numpy's astype(np.float32) does. NA / NaN stay NaN.
.as_float32 <- function(x) {
  v <- as.numeric(x)
  readBin(writeBin(v, raw(), size = 4L), what = "numeric", size = 4L, n = length(v))
}

# Feature extraction is implemented in R/staging_features.R.
