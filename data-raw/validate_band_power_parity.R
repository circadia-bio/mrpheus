# data-raw/validate_band_power_parity.R
#
# Compare mrpheus compute_band_power() against YASA's Python bandpower().
#
# Workflow
# --------
# 1. Generate the Python reference (once):
#
#      source /tmp/yasa_env/bin/activate
#      python3 data-raw/generate_bandpower_reference.py \
#          --edf ~/mne_data/physionet-sleep-data/SC4001E0-PSG.edf \
#          --channel "EEG Fpz-Cz" \
#          --out data-raw/yasa_reference_bandpower.csv
#
# 2. Set EDF_PATH and REF_CSV below and source this script.
#
# What to expect
# --------------
# Two known differences from YASA:
#   (a) Integration: YASA uses scipy.integrate.simpson, we use pracma::trapz.
#       Difference is typically < 0.01 % at nfft=1024.
#   (b) Band upper bound: YASA uses freq <= fmax (inclusive); we use freq < fmax
#       (exclusive). This means one frequency bin difference at each band edge.
#       Expected MAE: < 1 % for most bands; may be slightly higher for narrow
#       bands (sigma: 12-16 Hz).
#
# The validation script flags any band with mean absolute relative error (MARE)
# above MARE_THRESHOLD. Adjust as needed.

library(mrpheus)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

# ── Configuration ──────────────────────────────────────────────────────────────
EDF_PATH      <- "~/mne_data/physionet-sleep-data/SC4001E0-PSG.edf"
REF_CSV       <- "data-raw/yasa_reference_bandpower.csv"
EEG_CHANNEL   <- "EEG Fpz-Cz"
MARE_THRESHOLD <- 0.02   # flag bands with mean absolute relative error > 2 %

# ── Load EDF and compute mrpheus band power ────────────────────────────────────
message("Loading EDF with mrpheus...")
edf <- read_edf(EDF_PATH)
psg <- prepare_psg(edf, eeg_pattern = "EEG|Fpz|Pz")

message("Computing band power (absolute)...")
bp_abs <- compute_band_power(
  psg,
  channels  = EEG_CHANNEL,
  relative  = FALSE,
  win_sec   = 4,
  noverlap  = 200L,
  nfft      = 1024L
)

message("Computing band power (relative)...")
bp_rel <- compute_band_power(
  psg,
  channels  = EEG_CHANNEL,
  relative  = TRUE,
  win_sec   = 4,
  noverlap  = 200L,
  nfft      = 1024L
)

# Reshape to long format
bands <- c("delta", "theta", "alpha", "sigma", "beta", "gamma")

r_long_abs <- bp_abs |>
  select(epoch, channel, all_of(bands)) |>
  pivot_longer(all_of(bands), names_to = "band", values_to = "power_r")

r_long_rel <- bp_rel |>
  select(epoch, channel, all_of(bands)) |>
  pivot_longer(all_of(bands), names_to = "band", values_to = "relative_power_r")

r_long <- left_join(r_long_abs, r_long_rel, by = c("epoch", "channel", "band"))

# ── Load Python reference ──────────────────────────────────────────────────────
message("Loading Python reference...")
py_ref <- read_csv(REF_CSV, show_col_types = FALSE) |>
  rename(power_py          = power,
         relative_power_py = relative_power)

# ── Join and compare ───────────────────────────────────────────────────────────
# Limit to epochs present in both (Python may stop earlier if it hits EOF)
joined <- inner_join(
  r_long,
  py_ref |> select(epoch, band, power_py, relative_power_py),
  by = c("epoch", "band")
)

n_matched <- n_distinct(joined$epoch)
message(sprintf("Matched %d epochs across both implementations.", n_matched))

# ── Per-band summary ──────────────────────────────────────────────────────────
summary_tbl <- joined |>
  group_by(band) |>
  summarise(
    n             = n(),
    # Absolute power
    mae_abs       = mean(abs(power_r - power_py)),
    mare_abs      = mean(abs(power_r - power_py) / (abs(power_py) + 1e-20)),
    max_ae_abs    = max(abs(power_r - power_py)),
    # Relative power
    mae_rel       = mean(abs(relative_power_r - relative_power_py)),
    mare_rel      = mean(abs(relative_power_r - relative_power_py) /
                           (abs(relative_power_py) + 1e-20)),
    max_ae_rel    = max(abs(relative_power_r - relative_power_py)),
    .groups       = "drop"
  ) |>
  mutate(
    flagged = mare_rel > MARE_THRESHOLD
  ) |>
  arrange(band)

cat("\n── Band power parity report ──────────────────────────────────────────\n")
print(summary_tbl, n = Inf)

n_flagged <- sum(summary_tbl$flagged)
if (n_flagged == 0) {
  cat(sprintf("\n✓ All bands within %.0f %% MARE threshold.\n",
              MARE_THRESHOLD * 100))
} else {
  cat(sprintf("\n✖ %d band(s) exceed %.0f %% MARE threshold: %s\n",
              n_flagged,
              MARE_THRESHOLD * 100,
              paste(summary_tbl$band[summary_tbl$flagged], collapse = ", ")))
}

# ── Save report ────────────────────────────────────────────────────────────────
write_csv(summary_tbl, "data-raw/parity_report_bandpower.csv")
message("Saved: data-raw/parity_report_bandpower.csv")

# ── Scatter plots: R vs Python per band ───────────────────────────────────────
p <- joined |>
  ggplot(aes(x = relative_power_py, y = relative_power_r)) +
  geom_point(alpha = 0.15, size = 0.8) +
  geom_abline(slope = 1, intercept = 0, colour = "#B83E2C", linewidth = 0.6) +
  facet_wrap(~band, scales = "free") +
  labs(
    title    = "Band power parity: mrpheus vs YASA",
    subtitle = sprintf("SC4001E0-PSG.edf | %s | %d epochs",
                       EEG_CHANNEL, n_matched),
    x        = "YASA (Python)",
    y        = "mrpheus (R)"
  ) +
  theme_minimal(base_size = 11)

ggsave("data-raw/parity_bandpower.png", p,
       width = 10, height = 7, dpi = 150)
message("Saved: data-raw/parity_bandpower.png")
