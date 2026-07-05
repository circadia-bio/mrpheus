# dev/inspect_staging_model.R
#
# Inspects the bundled YASA LightGBM staging model and maps its 149-feature
# spec to the R implementation in .extract_staging_features().
#
# Run this interactively to verify the model loads correctly and to see
# exactly which features still need implementing.
#
# Reference: Vallat & Walker (2021), eLife 10:e70092

library(lightgbm)

model_path <- system.file("models", "yasa_staging.txt", package = "mrpheus")
# Or directly during development:
# model_path <- "inst/models/yasa_staging.txt"

model <- lgb.load(model_path)

# ── Model metadata ─────────────────────────────────────────────────────────────
cat("Classes     :", model$num_class(), "\n")        # expected: 5
cat("Features    :", model$num_feature(), "\n")      # expected: 149
cat("Class order : W  N1  N2  N3  REM\n")           # index 0–4

feat <- model$feature_name()
cat("\nAll feature names:\n")
print(feat)

# ── Feature groups ─────────────────────────────────────────────────────────────
#
# The 149 features break into four blocks:
#
#   EEG   (indices  0– 62):  21 base features × 3 variants = 63
#   time  (indices 63– 64):  2 features (hour, normalised position in night)
#   EOG   (indices 65–115): 17 base features × 3 variants = 51
#   EMG   (indices 116–148): 11 base features × 3 variants = 33
#
# Variant suffixes:
#   (none)        – raw value for that epoch
#   _c7min_norm   – z-scored relative to a ±7-min centred rolling window
#                   (zoo::rollapply over ±14 epochs)
#   _p2min_norm   – z-scored relative to the preceding 2-min window
#                   (zoo::rollapply over previous 4 epochs)

eeg_feats  <- feat[startsWith(feat, "eeg_")]
time_feats <- feat[startsWith(feat, "time_")]
eog_feats  <- feat[startsWith(feat, "eog_")]
emg_feats  <- feat[startsWith(feat, "emg_")]

cat("\nEEG features  (", length(eeg_feats),  "):", eeg_feats,  "\n")
cat("Time features (", length(time_feats), "):", time_feats, "\n")
cat("EOG features  (", length(eog_feats),  "):", eog_feats,  "\n")
cat("EMG features  (", length(emg_feats),  "):", emg_feats,  "\n")

stopifnot(length(feat) == 149L)
stopifnot(length(eeg_feats)  == 63L)
stopifnot(length(time_feats) ==  2L)
stopifnot(length(eog_feats)  == 51L)
stopifnot(length(emg_feats)  == 33L)

# ── Base features by signal type ───────────────────────────────────────────────
#
# Strip the channel prefix and variant suffix to recover the unique base names.
strip <- function(x, prefix) {
  x <- sub(paste0("^", prefix, "_"), "", x)
  x[!grepl("_c7min_norm$|_p2min_norm$", x)]
}

eeg_base  <- strip(eeg_feats,  "eeg")
eog_base  <- strip(eog_feats,  "eog")
emg_base  <- strip(emg_feats,  "emg")

cat("\nEEG base features (", length(eeg_base), "):\n"); print(eeg_base)
cat("\nEOG base features (", length(eog_base), "):\n"); print(eog_base)
cat("\nEMG base features (", length(emg_base), "):\n"); print(emg_base)

# ── Implementation status ──────────────────────────────────────────────────────
#
# What compute_band_power() already provides (relative = TRUE):
#   theta, alpha, sigma, beta  — correct bands
#   delta (0.5–4 Hz unified)   — WRONG: YASA splits into sdelta + fdelta
#   total_power                — this is abspow (use relative = FALSE)
#
# EEG spectral bands needed (YASA definitions):
#   sdelta : 0.5–2 Hz   (slow delta)
#   fdelta : 2–4 Hz     (fast delta)
#   theta  : 4–8 Hz
#   alpha  : 8–13 Hz    (note: 8–13 not 8–12)
#   sigma  : 12–16 Hz
#   beta   : 16–30 Hz
#   abspow : 0.5–40 Hz  (absolute total power, not relative)
#
# Spectral ratio features (EEG only):
#   dt = sdelta / theta
#   ds = sdelta / sigma
#   db = sdelta / beta
#   at = alpha  / theta
#
# Statistical time-domain features (EEG + EOG; EEG + EOG + EMG for nonlinear):
#   std      — standard deviation
#   iqr      — interquartile range
#   skew     — skewness (e1071::skewness or moments::skewness)
#   kurt     — excess kurtosis
#   nzc      — number of zero crossings (sum(diff(sign(x)) != 0))
#
# Nonlinear features (EEG + EOG + EMG):
#   perm     — permutation entropy (DescTools or custom)
#   higuchi  — Higuchi fractal dimension (pracma::fd_higuchi or custom)
#   hmob     — Hjorth mobility = sd(diff(x)) / sd(x)
#   hcomp    — Hjorth complexity = hmob(diff(x)) / hmob(x)
#   petrosian — Petrosian fractal dimension
#                = log10(N) / (log10(N) + log10(N / (N + 0.4 * nzc)))
#
# Context normalisation (all features except time):
#   _c7min_norm : (x - roll_mean(x, 15, fill="extend")) /
#                   (roll_sd(x, 15, fill="extend") + eps)
#                 window = 15 epochs = ±7 epochs = ±3.5 min ≈ "7 min centred"
#                 Use zoo::rollapply(x, 15, mean/sd, fill=NA, align="center")
#   _p2min_norm : (x - roll_mean(x, 4, fill="extend")) /
#                   (roll_sd(x, 4, fill="extend") + eps)
#                 window = 4 epochs = 2 min preceding
#                 Use zoo::rollapply(x, 4, mean/sd, fill=NA, align="right")
#
# eps = 1e-10 (matches YASA, avoids division by zero)
#
# Time features:
#   time_hour : elapsed hours since first epoch  (seq(0, n-1) * 30 / 3600)
#   time_norm : time_hour / max(time_hour)        (0 → 1 across the recording)

cat("\n── TODO in .extract_staging_features() ──────────────────────────────────\n")
cat("[ ] Fix band defs: split delta -> sdelta (0.5-2) + fdelta (2-4)\n")
cat("[ ] Change alpha upper bound: 8-13 Hz (not 8-12)\n")
cat("[ ] abspow: call compute_band_power(relative=FALSE), sum all bands\n")
cat("[ ] Spectral ratios: dt, ds, db, at\n")
cat("[ ] Statistical features: std, iqr, skew, kurt, nzc\n")
cat("[ ] Nonlinear features: perm, higuchi, hmob, hcomp, petrosian\n")
cat("[ ] EOG feature set (no ratios, rest same as EEG)\n")
cat("[ ] EMG feature set (abspow + nonlinear only, no spectral bands)\n")
cat("[ ] Time features: time_hour, time_norm\n")
cat("[ ] Context normalisation: _c7min_norm (zoo, align=center, k=15)\n")
cat("[ ]                        _p2min_norm (zoo, align=right,  k=4)\n")
cat("[ ] Final column order must match feat vector exactly\n")
