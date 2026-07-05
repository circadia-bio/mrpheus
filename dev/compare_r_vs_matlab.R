# dev/compare_r_vs_matlab.R
#
# Visual comparison of detect_qrs() (R) vs pan_tompkin.m (MATLAB) outputs.
# Run from the mrpheus package root after generating MATLAB fixtures:
#   source("dev/compare_r_vs_matlab.R")

library(mrpheus)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gsignal)

FIXTURE_DIR <- "tests/testthat/fixtures"

# ---- Check fixtures exist --------------------------------------------------
needed  <- c("fixture_qrs_i_raw.csv", "fixture_qrs_amp_raw.csv")
missing <- needed[!file.exists(file.path(FIXTURE_DIR, needed))]
if (length(missing) > 0)
  stop("Fixtures not found. Run: matlab -batch \"run data-raw/generate_matlab_fixtures.m\"\n",
       "Missing: ", paste(missing, collapse = ", "))

# ---- Load data -------------------------------------------------------------
path  <- system.file("extdata", "example_physlog.log", package = "mrpheus")
rec   <- read_philips_physlog(path)
ecg   <- as.double(rec$C[, "v1raw"])
sfreq <- rec$HDR$sfreq

# R output
qrs_r <- detect_qrs(ecg, fs = sfreq)

# MATLAB output
matlab_i   <- as.integer(read.csv(file.path(FIXTURE_DIR, "fixture_qrs_i_raw.csv"),
                                   header = FALSE)[[1]])
matlab_amp <- as.double(read.csv(file.path(FIXTURE_DIR, "fixture_qrs_amp_raw.csv"),
                                  header = FALSE)[[1]])

# Bandpass-filtered signal for plotting
bp     <- gsignal::butter(3L, c(5, 15) * 2 / sfreq, type = "pass")
ecg_bp <- as.double(gsignal::filtfilt(bp, ecg))
t      <- seq_along(ecg_bp) / sfreq

cat("── Peak count ──────────────────────────────\n")
cat("R     :", length(qrs_r$qrs_i), "peaks\n")
cat("MATLAB:", length(matlab_i),    "peaks\n")

theme_qc <- function() {
  theme_minimal(base_size = 13) +
    theme(
      panel.grid        = element_blank(),
      axis.line         = element_line(colour = "black"),
      axis.line.x.top   = element_blank(),
      axis.line.y.right = element_blank(),
      legend.position   = "top"
    )
}

# ---- Plot 1: ECG with overlaid peaks ---------------------------------------
df_ecg <- data.frame(t = t, ecg = ecg_bp)

df_peaks <- bind_rows(
  data.frame(t = qrs_r$qrs_i / sfreq, amp = qrs_r$qrs_amp,  src = "R"),
  data.frame(t = matlab_i   / sfreq,  amp = matlab_amp,      src = "MATLAB")
)

p1 <- ggplot(df_ecg, aes(t, ecg)) +
  geom_line(linewidth = 0.3, colour = "grey50") +
  geom_point(data = df_peaks,
             aes(x = t, y = amp, colour = src, shape = src),
             size = 2.5, alpha = 0.85) +
  scale_colour_manual(values = c(R = "#014370", MATLAB = "#FC544A")) +
  scale_shape_manual(values  = c(R = 19,        MATLAB = 4)) +
  labs(title = "R-peak detection: R vs MATLAB",
       x = "Time (s)", y = "Amplitude (normalised)",
       colour = NULL, shape = NULL) +
  theme_qc()

print(p1)

# ---- Plot 2: RR intervals over time ----------------------------------------
df_rr <- bind_rows(
  data.frame(t = qrs_r$qrs_i[-1] / sfreq,
             rr = diff(qrs_r$qrs_i) / sfreq * 1000, src = "R"),
  data.frame(t = matlab_i[-1] / sfreq,
             rr = diff(matlab_i)   / sfreq * 1000,   src = "MATLAB")
)

p2 <- ggplot(df_rr, aes(t, rr, colour = src)) +
  geom_line(linewidth = 0.6, alpha = 0.8) +
  geom_point(size = 1.5, alpha = 0.7) +
  scale_colour_manual(values = c(R = "#014370", MATLAB = "#FC544A")) +
  labs(title = "RR intervals: R vs MATLAB",
       x = "Time (s)", y = "RR interval (ms)", colour = NULL) +
  theme_qc()

print(p2)

# ---- Plot 3: Sample-by-sample index difference -----------------------------
if (length(qrs_r$qrs_i) == length(matlab_i)) {

  df_diff <- data.frame(
    peak  = seq_along(qrs_r$qrs_i),
    delta = sort(qrs_r$qrs_i) - sort(matlab_i)
  )

  p3 <- ggplot(df_diff, aes(peak, delta)) +
    geom_hline(yintercept = 0,      linetype = "dashed", colour = "grey60") +
    geom_hline(yintercept = c(-5, 5), linetype = "dotted",
               colour = "#FC544A", alpha = 0.6) +
    geom_point(size = 1.8, colour = "#014370", alpha = 0.8) +
    geom_line(linewidth = 0.3, colour = "#014370", alpha = 0.5) +
    labs(title    = "R-peak index difference (R minus MATLAB)",
         subtitle = "Dotted lines = ±5 sample tolerance",
         x = "Peak number", y = "Difference (samples)") +
    theme_qc()

  print(p3)

  cat("\n── Index differences (samples) ─────────────\n")
  cat("Max :", max(abs(df_diff$delta)), "\n")
  cat("Mean:", round(mean(abs(df_diff$delta)), 2), "\n")
  cat("All within ±5:", all(abs(df_diff$delta) <= 5), "\n")

} else {
  cat("\nPeak counts differ (R:", length(qrs_r$qrs_i),
      "vs MATLAB:", length(matlab_i), ") — skipping difference plot.\n")
}
