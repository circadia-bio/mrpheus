# data-raw/create_example_physlog.R
#
# Generates inst/extdata/example_physlog.log — a synthetic Philips wBTU
# physiological log for use in the mri-physiology vignette.
#
# Run once from the package root:
#   source("data-raw/create_example_physlog.R")
#
# Signal characteristics:
#   - Duration  : 20 s at 496 Hz (9920 samples)
#   - ECG       : neonatal HR ~150 bpm, Gaussian QRS complexes
#   - PPU       : cardiac pulsation with ~150 ms peripheral delay
#   - RESP      : ~40 breaths/min (neonatal)
#   - Accel     : near-zero (subject still)
#   - Markers   : ScannerStart at 1 s, ScannerStop near end, VcgOnset at
#                 each R-peak

set.seed(42)

sfreq    <- 496L
duration <- 20L
n        <- sfreq * duration

t <- seq(0, duration - 1 / sfreq, by = 1 / sfreq)

# ---- ECG (v1raw) ------------------------------------------------------------
# Neonatal HR ~150 bpm -> RR interval = 198 samples
rr_samp <- round(60 / 150 * sfreq)
r_peaks <- seq(rr_samp, n, by = rr_samp)

ecg_raw <- numeric(n)
for (r in r_peaks) {
  win <- max(1L, r - 120L):min(n, r + 180L)
  rel <- win - r
  ecg_raw[win] <- ecg_raw[win] +
    -250 * exp(-(rel + 12L)^2 / (2 * 4^2))  +   # Q wave
    1800 * exp(-rel^2          / (2 * 3^2))  +   # R wave
    -300 * exp(-(rel - 10L)^2  / (2 * 4^2)) +   # S wave
     400 * exp(-(rel - 45L)^2  / (2 * 18^2))    # T wave
}
ecg_raw <- ecg_raw + 50 * sin(2 * pi * 0.1 * t) + rnorm(n, 0, 30)

v1raw <- as.integer(round(ecg_raw))
v2raw <- as.integer(round(ecg_raw * 0.6  + rnorm(n, 0, 20)))
v1    <- as.integer(round(ecg_raw * 0.02))
v2    <- as.integer(round(ecg_raw * 0.012))

# ---- PPU --------------------------------------------------------------------
dppu    <- round(0.15 * sfreq)  # ~150 ms peripheral delay
ppu_raw <- numeric(n)
for (r in r_peaks) {
  rd  <- r + dppu
  win <- max(1L, rd - 80L):min(n, rd + 120L)
  ppu_raw[win] <- ppu_raw[win] + 600 * exp(-(win - rd)^2 / (2 * 28^2))
}
ppu <- as.integer(round(ppu_raw - 300 + rnorm(n, 0, 15)))

# ---- RESP -------------------------------------------------------------------
resp <- as.integer(round(500 * sin(2 * pi * (40 / 60) * t) + rnorm(n, 0, 20)))

# ---- Accelerometer ----------------------------------------------------------
gx <- as.integer(round(rnorm(n, 0, 2)))
gy <- as.integer(round(rnorm(n, 0, 2)))
gz <- as.integer(round(rnorm(n, 0, 2)))

# ---- Markers ----------------------------------------------------------------
mark <- integer(n)
mark[sfreq]               <- bitwOr(mark[sfreq],           16L)  # ScannerStart
mark[n - sfreq + 42L]     <- bitwOr(mark[n - sfreq + 42L], 32L)  # ScannerStop
for (r in r_peaks)
  if (r >= 1L && r <= n) mark[r] <- bitwOr(mark[r], 1L)          # VcgOnset

mark_hex <- sprintf("%04X", mark)

# ---- Write ------------------------------------------------------------------
out_dir <- file.path("inst", "extdata")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_path <- file.path(out_dir, "example_physlog.log")

con <- file(out_path, "w")
writeLines(c(
  "## Northumbria University, Release r32 (SWID 119)",
  "## 01-01-2024 10:00:00",
  "## -1071 377 -592 257 1618 -697 605 -251 0",
  "## Dockable table = FALSE",
  "# v1raw v2raw  v1 v2  ppu resp  gx gy gz mark"
), con)

for (i in seq_len(n)) {
  writeLines(
    paste(v1raw[i], v2raw[i], v1[i], v2[i], ppu[i], resp[i],
          gx[i], gy[i], gz[i], mark_hex[i]),
    con
  )
}
close(con)

message(
  "Written: ", out_path,
  " (", n, " samples / ", duration, " s / ",
  length(r_peaks), " R-peaks)"
)
