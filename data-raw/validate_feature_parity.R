# data-raw/validate_feature_parity.R
#
# Compare mrpheus feature extraction against YASA Python reference output.
#
# Workflow
# --------
# 1. Generate the Python reference (once, from the same EDF):
#
#      source /tmp/yasa_env/bin/activate
#      python3 data-raw/generate_yasa_reference.py \
#          --edf /path/to/recording.edf \
#          --eeg "EEG Fpz-Cz" \
#          --eog "EOG horizontal" \
#          --emg "EMG submental" \
#          --out data-raw/yasa_reference_features.csv
#
# 2. Set the two paths below and source this script in R.
#
# The script reports mean absolute error (MAE) and max absolute error per
# feature, flags any columns exceeding the tolerance threshold, and writes a
# tidy summary CSV to data-raw/parity_report.csv.

library(mrpheus)
library(dplyr)
library(tidyr)
library(ggplot2)

# ── Configuration ─────────────────────────────────────────────────────────────
EDF_PATH     <- "/path/to/recording.edf"    # <-- set to your EDF
REF_CSV      <- "data-raw/yasa_reference_features.csv"
EEG_CHANNEL  <- "EEG Fpz-Cz"               # <-- match what you passed to Python
EOG_CHANNEL  <- "EOG horizontal"            # <-- or NA_character_ if omitted
EMG_CHANNEL  <- "EMG submental"             # <-- or NA_character_ if omitted

# Tolerance: raw features should agree to this level (relative powers are
# unitless 0-1; absolute power in µV²/Hz may have larger absolute error).
MAE_THRESHOLD <- 0.01   # flag features with MAE above this

# Edge epochs to exclude from comparison (rolling normalisation differs at
# boundaries; both use partial windows but may handle fill differently).
EDGE_EPOCHS <- 15L   # skip first and last N epochs

# ── Load EDF and extract mrpheus features ────────────────────────────────────
message("Loading EDF with mrpheus...")
edf <- read_edf(EDF_PATH)
psg <- prepare_psg(edf)

message("Extracting features...")
r_feats <- mrpheus:::.extract_staging_features(
  psg, EEG_CHANNEL, EOG_CHANNEL, EMG_CHANNEL
)

# Drop the epoch index column for comparison
r_mat <- r_feats |> select(-epoch) |> as.matrix()

# ── Load Python reference ─────────────────────────────────────────────────────
message("Loading YASA reference features...")
py_feats <- read.csv(REF_CSV, check.names = FALSE)
py_mat   <- as.matrix(py_feats)

# ── Structural checks ─────────────────────────────────────────────────────────
message("\n── Structural checks ──────────────────────────────────────────────────")

cat("mrpheus epochs  :", nrow(r_mat), "\n")
cat("YASA    epochs  :", nrow(py_mat), "\n")
cat("mrpheus features:", ncol(r_mat), "\n")
cat("YASA    features:", ncol(py_mat), "\n")

if (ncol(r_mat) != ncol(py_mat)) {
  stop("Feature count mismatch: mrpheus=", ncol(r_mat), " vs YASA=", ncol(py_mat))
}

# Column name comparison
r_names  <- colnames(r_mat)
py_names <- colnames(py_mat)

if (!identical(r_names, py_names)) {
  missing_from_r  <- setdiff(py_names, r_names)
  extra_in_r      <- setdiff(r_names, py_names)
  order_mismatch  <- r_names[r_names != py_names]

  message("\n⚠  Column name issues detected:")
  if (length(missing_from_r) > 0)
    message("  Missing from mrpheus : ", paste(missing_from_r, collapse = ", "))
  if (length(extra_in_r) > 0)
    message("  Extra in mrpheus     : ", paste(extra_in_r, collapse = ", "))
  if (length(order_mismatch) > 0)
    message("  Order mismatch at    : ", paste(order_mismatch[1:min(5, length(order_mismatch))], collapse = ", "))

  # Reorder Python columns to match mrpheus for numeric comparison
  common <- intersect(r_names, py_names)
  r_mat  <- r_mat[,  common, drop = FALSE]
  py_mat <- py_mat[, common, drop = FALSE]
  message("  Proceeding with ", length(common), " common columns.\n")
} else {
  message("✔  Column names match exactly (", ncol(r_mat), " features)\n")
}

# Align epoch counts (use minimum)
n_epochs <- min(nrow(r_mat), nrow(py_mat))
r_mat  <- r_mat[seq_len(n_epochs),  , drop = FALSE]
py_mat <- py_mat[seq_len(n_epochs), , drop = FALSE]

# ── Numerical comparison (excluding edge epochs) ──────────────────────────────
inner <- seq(EDGE_EPOCHS + 1L, n_epochs - EDGE_EPOCHS)
r_inner  <- r_mat[inner,  , drop = FALSE]
py_inner <- py_mat[inner, , drop = FALSE]

diff_mat <- abs(r_inner - py_inner)

report <- tibble::tibble(
  feature = colnames(diff_mat),
  mae     = colMeans(diff_mat, na.rm = TRUE),
  max_ae  = apply(diff_mat, 2, max, na.rm = TRUE),
  flagged = mae > MAE_THRESHOLD
)

# ── Print summary ─────────────────────────────────────────────────────────────
message("── Parity report (inner ", length(inner), " epochs, edge ±", EDGE_EPOCHS, " excluded) ──────")

n_flagged <- sum(report$flagged)
if (n_flagged == 0) {
  message("✔  All ", nrow(report), " features within MAE threshold (", MAE_THRESHOLD, ")\n")
} else {
  message("⚠  ", n_flagged, " / ", nrow(report), " features exceed MAE threshold:\n")
  report |>
    filter(flagged) |>
    arrange(desc(mae)) |>
    print(n = Inf)
}

message("\nTop 10 features by MAE:")
report |>
  arrange(desc(mae)) |>
  slice_head(n = 10) |>
  mutate(across(c(mae, max_ae), \(x) round(x, 6))) |>
  print(n = Inf)

# ── Visualise ─────────────────────────────────────────────────────────────────
# For flagged features, plot R vs Python across epochs side by side
if (n_flagged > 0) {
  flagged_feats <- report$feature[report$flagged]
  n_plot        <- min(6L, length(flagged_feats))

  plot_data <- purrr::map_dfr(flagged_feats[seq_len(n_plot)], function(f) {
    tibble::tibble(
      epoch   = inner,
      feature = f,
      r_value = r_inner[, f],
      py_value = py_inner[, f]
    )
  }) |>
    tidyr::pivot_longer(c(r_value, py_value),
                        names_to  = "source",
                        values_to = "value") |>
    mutate(source = if_else(source == "r_value", "mrpheus", "YASA"))

  p <- ggplot(plot_data, aes(x = epoch, y = value, colour = source)) +
    geom_line(linewidth = 0.4, alpha = 0.8) +
    facet_wrap(~ feature, scales = "free_y", ncol = 2) +
    scale_colour_manual(values = c(mrpheus = "#014370", YASA = "#FC544A")) +
    labs(title = "Flagged features: mrpheus vs YASA",
         x = "Epoch", y = "Feature value", colour = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom",
          panel.grid.minor = element_blank())

  ggsave("data-raw/parity_flagged_features.png", p,
         width = 10, height = 2.5 * ceiling(n_plot / 2), dpi = 150)
  message("Plot saved: data-raw/parity_flagged_features.png")
}

# ── Save report ───────────────────────────────────────────────────────────────
readr::write_csv(report, "data-raw/parity_report.csv")
message("Report saved: data-raw/parity_report.csv")
