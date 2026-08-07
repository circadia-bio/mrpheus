## mrpheus 0.1.4.9000  (development)

### 🚀 Performance fixes in `detect_qrs()`

`detect_qrs()` was hanging indefinitely on real recordings (unbounded on a
468,700-sample / ~15.7 min recording at ~500 Hz). Two separate O(n²)
bottlenecks were found and fixed:

* **Peak selection** (`.find_peaks_min_dist()`): was checking every candidate
  peak against every other accepted peak on each iteration — fine on toy
  data, unusable once candidate counts reach the hundreds of thousands.
  Replaced with a binary search (`findInterval()`) against a pre-allocated,
  sorted buffer of accepted peak locations. Note this is O(n_candidates log k
  + k²) where k is the final peak count (bounded by physiological heart
  rate, not by recording length), rather than a universal O(n log n) — fine
  here, but worth remembering if this helper is ever reused somewhere the
  output peak count could scale with input size.
* **Moving-average integration** (Stage 4): used `stats::convolve()` for a
  plain boxcar filter. R's FFT-based convolution can degrade badly depending
  on how the padded length factors for `stats::nextn()`, and this was the
  actual bottleneck once the peak-selection fix was in place. Replaced with a
  cumulative-sum boxcar implementation, mathematically equivalent to the
  FFT-based full convolution (verified numerically identical, differences on
  the order of 1e-15) and unconditionally O(n).

Combined effect: `detect_qrs()` on the 468,700-sample recording above now
completes in 0.97 seconds (1855 R-peaks detected, mean HR 118 bpm).

### ✨ New function: `compute_hrv_freq()`

* `compute_hrv_freq()` — frequency-domain HRV (HF power) locked to each
  infant's own respiratory rate, rather than the fixed adult-derived HF band
  (0.15–0.4 Hz) that doesn't apply to neonatal breathing rates (~40–70
  breaths/min at rest). Detects the dominant respiratory frequency directly
  from the physlog's `resp` channel, resamples the RR-interval series via
  cubic spline, and integrates Welch PSD power in a band centered on that
  detected frequency. Verified against synthetic data: exact frequency
  recovery and correct amplitude-scaling of integrated power.
* Fills the gap left by `compute_hrv_sleep()`, which remains a stub expecting
  staged PSG input (wrong shape for physlog data) with its RR-to-HRV step
  never wired in — `compute_hrv_freq()` is a separate, physlog-shaped path
  rather than a completion of that stub.

### 🧪 Tests

* `test-find-peaks-min-dist.R`: `.find_peaks_min_dist()` against a
  brute-force reference on random and periodic synthetic signals, the
  minimum-distance constraint, tallest-peak-wins tie-breaking (including the
  all-tied flat-signal edge case), and a timing regression guard on a
  large, realistically-smoothed candidate set.
* `test-boxcar-integration.R`: the cumulative-sum boxcar formula against
  `stats::convolve(type = "open")` across several window sizes and signal
  lengths, plus an end-to-end `detect_qrs()` timing smoke test on a
  synthetic two-minute ECG-like signal.
* `test-compute_hrv_freq.R`: synthetic QRS/respiration generators injecting
  a known respiratory frequency via RSA-modulated RR intervals — frequency
  recovery at two different frequencies, HF power scaling with RSA
  amplitude, band-bounds correctness, the too-few-beats `NA` path, and
  implausible-RR-interval filtering.

All fixtures are fully synthetic — no participant-derived data.

## mrpheus 0.1.4  (2026-07)

### 🚀 Rcpp hot-path optimisation (staging pipeline + event detection)

Complete C++ hot-path port across the staging and event detection pipelines,
replacing the remaining R-level bottlenecks. `zoo` removed from Imports.

**Staging pipeline — per-epoch functions** (`src/staging_features.cpp`):

- `nzc_cpp`: zero-crossing count — single O(N) pass replacing `sum(diff(sign(x)) != 0)`
- `petrosian_fd_cpp`: local extrema count + log formula — replaces `diff` + logical sum
- `hjorth_cpp`: mobility and complexity — single O(N) pass accumulating all three variances
  (x, diff(x), diff(diff(x))) without materialising intermediate vectors
- `stat_features_cpp`: std, IQR (type-7 quantile), skewness, excess kurtosis — replaces
  `sd()`, `IQR()` (which sorts), and two `mean()` calls per epoch
- `rowmedian_cpp`: row-wise median of the Welch periodogram matrix — eliminates
  ~750 000 R `median()` dispatch calls per recording (251 freq bins x ~3 000 epochs)
- `roll_right_mean_cpp`: right-aligned partial rolling mean (k=4), replacing
  `zoo::rollapply` for the `_p2min_norm` normalisation step (49 calls/recording)
- `robust_scale_cpp`: `(x - median) / (q95 - q5)` with type-7 quantile and
  NA-awareness, replacing `median()` + `quantile()` for 49 normalisation calls

**Event detection** (`src/event_detection.cpp`, new file):

- `roll_rms_cpp`: centered rolling RMS replacing `zoo::rollapply` in
  `compute_spindles()` — processes the full concatenated epoch signal per channel
- `detect_so_candidates_cpp`: zero-crossing scan + duration/amplitude filtering
  for slow oscillation detection, replacing the R `which()` + `lapply` loop
  in `compute_slow_oscillations()`

### 🔧 Other changes

- Fixed a latent bug where the pure-R `vapply` definition of `.roll_triang_mean`
  in the normalisation section was shadowing `roll_triang_mean_cpp` on every
  `load_all()`. Now a single Rcpp-backed definition exists.
- Removed `zoo` from Imports — no remaining `zoo::` calls anywhere in the package.

### 🧪 Tests

All 371 tests pass. `0 errors | 0 warnings | 0 notes`.

## mrpheus 0.1.3  (2026-07)

### 🩺 PSG preprocessing pipeline

Full preprocessing pipeline from EDF to analysis-ready epochs, matching the
companion Python/MNE notebook workflow.

* `preprocess_psg()` — continuous-signal preprocessing pipeline: channel
  renaming, DC removal, powerline notch filter (auto-detected 50/60 Hz +
  harmonics), and channel-type-specific bandpass filtering (EEG/EOG/EMG/ECG).
  Operates on the full continuous signal before re-epoching to avoid filter
  discontinuities at epoch boundaries.
* `remove_dc()`, `detect_powerline()`, `notch_filter()`, `bandpass_filter()` —
  exported signal-level functions; each operates on a plain numeric vector
  independently of the PSG pipeline.
