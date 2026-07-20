# Changelog

## mrpheus 0.1.4

#### Rcpp hot-path optimisation (staging pipeline + event detection)

Complete C++ hot-path port across the staging and event detection
pipelines, replacing the remaining R-level bottlenecks. `zoo` removed
from Imports.

**Staging pipeline — per-epoch functions** (`src/staging_features.cpp`):

- `nzc_cpp`: zero-crossing count — single O(N) pass replacing
  `sum(diff(sign(x)) != 0)`
- `petrosian_fd_cpp`: local extrema count + log formula — replaces
  `diff` + logical sum
- `hjorth_cpp`: mobility and complexity — single O(N) pass accumulating
  all three variances (x, diff(x), diff(diff(x))) without materialising
  intermediate vectors
- `stat_features_cpp`: std, IQR (type-7 quantile), skewness, excess
  kurtosis — replaces [`sd()`](https://rdrr.io/r/stats/sd.html),
  [`IQR()`](https://rdrr.io/r/stats/IQR.html) (which sorts), and two
  [`mean()`](https://rdrr.io/r/base/mean.html) calls per epoch
- `rowmedian_cpp`: row-wise median of the Welch periodogram matrix —
  eliminates ~750 000 R
  [`median()`](https://rdrr.io/r/stats/median.html) dispatch calls per
  recording (251 freq bins x ~3 000 epochs)
- `roll_right_mean_cpp`: right-aligned partial rolling mean (k=4),
  replacing `zoo::rollapply` for the `_p2min_norm` normalisation step
  (49 calls/recording)
- `robust_scale_cpp`: `(x - median) / (q95 - q5)` with type-7 quantile
  and NA-awareness, replacing
  [`median()`](https://rdrr.io/r/stats/median.html) +
  [`quantile()`](https://rdrr.io/r/stats/quantile.html) for 49
  normalisation calls

**Event detection** (`src/event_detection.cpp`, new file):

- `roll_rms_cpp`: centered rolling RMS replacing `zoo::rollapply` in
  [`compute_spindles()`](https://mrpheus.circadia-lab.uk/reference/compute_spindles.md)
  — processes the full concatenated epoch signal per channel
- `detect_so_candidates_cpp`: zero-crossing scan + duration/amplitude
  filtering for slow oscillation detection, replacing the R
  [`which()`](https://rdrr.io/r/base/which.html) + `lapply` loop in
  [`compute_slow_oscillations()`](https://mrpheus.circadia-lab.uk/reference/compute_slow_oscillations.md)

**Other changes:**

- Fixed a latent bug where the pure-R `vapply` definition of
  `.roll_triang_mean` in the normalisation section was shadowing
  `roll_triang_mean_cpp` on every `load_all()`. Now a single Rcpp-backed
  definition exists.
- Removed `zoo` from Imports — no remaining `zoo::` calls anywhere in
  the package.

All 371 tests pass. `0 errors | 0 warnings | 0 notes`.

## mrpheus 0.1.3

#### PSG preprocessing pipeline

Full preprocessing pipeline from EDF to analysis-ready epochs, matching
the companion Python/MNE notebook workflow.

- [`preprocess_psg()`](https://mrpheus.circadia-lab.uk/reference/preprocess_psg.md)
  — continuous-signal preprocessing pipeline: channel renaming, DC
  removal, powerline notch filter (auto-detected 50/60 Hz + harmonics),
  and channel-type-specific bandpass filtering (EEG/EOG/EMG/ECG).
  Operates on the full continuous signal before re-epoching to avoid
  filter discontinuities at epoch boundaries.
- [`remove_dc()`](https://mrpheus.circadia-lab.uk/reference/remove_dc.md),
  [`detect_powerline()`](https://mrpheus.circadia-lab.uk/reference/detect_powerline.md),
  [`notch_filter()`](https://mrpheus.circadia-lab.uk/reference/notch_filter.md),
  [`bandpass_filter()`](https://mrpheus.circadia-lab.uk/reference/bandpass_filter.md)
  — exported signal-level functions; each operates on a plain numeric
  vector independently of the PSG pipeline.
