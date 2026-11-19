# Changelog

## mrpheus 0.1.4.9000 (development)

### 🛏️ `stage_epochs()`: now reproduces YASA (results change)

[`stage_epochs()`](https://mrpheus.circadia-lab.uk/reference/stage_epochs.md)
is checked against `yasa.SleepStaging` (YASA 0.8.1, scipy 1.18.1)
feature by feature and probability by probability. Three errors in the
feature code, a missing model selection and a missing precision step are
fixed; on the public night used to check this, the stage is now the same
as YASA’s in every epoch.

- **Channel-specific classifiers.** YASA picks one of eight classifiers
  from the channels given (EEG, EEG + EOG, EEG + EMG, EEG + EOG + EMG,
  each with or without demographics). mrpheus shipped only the EEG +
  EOG + EMG model and filled missing channels with NA. With EOG absent
  that gave 50% accuracy against manual scoring on the test night (YASA:
  83%); with EMG absent it differed from YASA’s EEG + EOG model in 91 of
  780 epochs. All eight are now bundled as
  `inst/models/yasa_<channels>[_demo].txt.gz` (about 3.3 MB each) and
  chosen automatically. New arguments `age` and `male` select the
  `+demo` classifiers; `model_path` now defaults to `NULL` (automatic)
  and accepts `.gz` files. The result carries the classifier as
  `attr(, "classifier")`. Run `data-raw/fetch_yasa_models.py` to create
  the model files; the old `inst/models/yasa_staging.txt` is
  byte-identical to `yasa_eeg_eog_emg` and can be deleted.
- **Single precision.** YASA rounds all features to float32 before
  predicting; mrpheus now does too.
- **Simpson integration for an even number of points.**
  `.scipy_simpson()` used a 3/8-rule variant; scipy (\>= 1.11, checked
  on 1.14.1 and 1.18.1) uses Cartwright’s correction on the last
  interval. This changed the 0.4-1 Hz and 1-4 Hz band powers (even point
  counts at 0.2 Hz resolution) and the delta ratios built from them, by
  up to half a standard deviation, and also affects
  [`compute_band_power()`](https://mrpheus.circadia-lab.uk/reference/compute_band_power.md)
  and
  [`compute_temporal_bandpower()`](https://mrpheus.circadia-lab.uk/reference/compute_temporal_bandpower.md)
  wherever a band holds an even number of frequency bins. The old tests
  used linear data, which every Simpson variant integrates exactly, so
  they could not see it.
- **Skewness and kurtosis** standardised with the sample (n - 1)
  standard deviation; YASA uses population central moments. The
  normalised variants cancel the difference, the raw ones (which the
  classifiers also use) did not. A flat epoch now gives `NaN` skew and
  kurtosis, as in YASA.
- **Zero crossings** now follow `antropy.num_zerocross` (`np.signbit`):
  an exact zero is no longer counted as two crossings.

*Measured agreement.* On Sleep-EDF `SC4001E0` (6.5 h, 780 epochs), 100
Hz and 256 Hz inputs, all 16 classifier configurations give the same
stage in every epoch as YASA, with probabilities equal to about 1e-16 in
99.9% of epochs. The exception is the recording’s first and last sample:
the zero-phase FIR filter’s output there is zero by construction and its
sign is floating-point noise (+6e-15 in MNE, -4e-14 here), so one
zero-crossing count can differ by one. Because `robust_scale` uses
whole-column quantiles this shifts every epoch’s normalised value
slightly; the largest probability difference seen was 1e-3, in one epoch
of 12,480. Taking only the zero-crossing features from YASA makes the
probabilities exactly identical (difference 0.0).
`tests/testthat/test-yasa-parity.R` checks all eight classifiers against
a stored YASA reference. The earlier “validated” parity script compared
at MAE 0.01 and skipped 15 edge epochs; that was not parity to float
precision.

*Caveats.* Agreement is with YASA 0.8.1 on scipy \>= 1.15; older scipy
versions integrate even-length bands differently and give different YASA
features. Preprocessing with
[`preprocess_psg()`](https://mrpheus.circadia-lab.uk/reference/preprocess_psg.md)
before staging is not what YASA does and changes the result (84 of 780
epochs on the test night, mostly through the EOG band-pass); stage from
unpreprocessed signals to reproduce YASA.

### ⚠️ `compute_slow_oscillations()`: fixes and new defaults (results change)

- **Polarity bug fixed.** `detect_so_candidates_cpp()` paired zero
  crossings in steps of two starting from the first crossing in the
  epoch, whatever its direction. In an epoch whose filtered signal
  started negative every pair was half a cycle out of phase and the
  epoch returned no slow oscillations, so roughly half the epochs could
  be lost. The walk now starts at a crossing after which the signal is
  negative. Existing tests passed because their signals started at
  exactly 0.
- **Filter fixed at higher sample rates.** The single 8th-order
  band-pass was numerically singular at 256 Hz
  ([`gsignal::filter_zi()`](https://rdrr.io/pkg/gsignal/man/filter_zi.html)
  failed with “system is computationally singular”). It is now a
  zero-phase high-pass / low-pass cascade, designed once rather than per
  epoch.
- **New defaults, taken from `yasa.sw_detect()`:**
  `freq_so = c(0.3, 1.5)` (was `c(0.5, 2)`),
  `duration_neg_s = c(0.3, 1.5)` (was `c(0.1, 1.5)`),
  `amp_ptp_threshold_uv = c(75, 350)` (was `c(75, 500)`). New arguments
  `amp_neg_uv = c(40, 200)` and `amp_pos_uv = c(10, 150)` (trough and
  peak limits; `NULL` switches them off) and `transition_hz = 0.2`,
  which places the filter’s -6 dB points at 0.2 and 1.6 Hz as YASA’s FIR
  does. The old behaviour can be requested with `freq_so = c(0.5, 2)`,
  `duration_neg_s = c(0.1, 1.5)`, `amp_ptp_threshold_uv = c(75, 500)`,
  `amp_neg_uv = NULL`, `amp_pos_uv = NULL`, `transition_hz = 0`.
- **Not identical to YASA.** The pairing is still zero-crossing based,
  and epochs are filtered one at a time. In a Python comparison on one
  Sleep-EDF night (SC4001E0, Fpz-Cz, N2 + N3), a mirror of the old
  defaults gave 8.7 events per NREM minute against 5.9 for YASA, with
  53% of its events matching a YASA event (negative peaks within 0.5 s).
  A mirror of the new defaults gave 4.7 per minute, with 90% matching a
  YASA event and 72% of YASA’s events found. This is one night, one
  channel and a Python mirror, not a validation.

### 🚀 Performance fixes in `detect_qrs()`

[`detect_qrs()`](https://mrpheus.circadia-lab.uk/reference/detect_qrs.md)
was hanging indefinitely on real recordings (unbounded on a
468,700-sample / ~15.7 min recording at ~500 Hz). Two separate O(n²)
bottlenecks were found and fixed:

- **Peak selection** (`.find_peaks_min_dist()`): was checking every
  candidate peak against every other accepted peak on each iteration —
  fine on toy data, unusable once candidate counts reach the hundreds of
  thousands. Replaced with a binary search
  ([`findInterval()`](https://rdrr.io/r/base/findInterval.html)) against
  a pre-allocated, sorted buffer of accepted peak locations. Note this
  is O(n_candidates log k
  - k²) where k is the final peak count (bounded by physiological heart
    rate, not by recording length), rather than a universal O(n log n) —
    fine here, but worth remembering if this helper is ever reused
    somewhere the output peak count could scale with input size.
- **Moving-average integration** (Stage 4): used
  [`stats::convolve()`](https://rdrr.io/r/stats/convolve.html) for a
  plain boxcar filter. R’s FFT-based convolution can degrade badly
  depending on how the padded length factors for
  [`stats::nextn()`](https://rdrr.io/r/stats/nextn.html), and this was
  the actual bottleneck once the peak-selection fix was in place.
  Replaced with a cumulative-sum boxcar implementation, mathematically
  equivalent to the FFT-based full convolution (verified numerically
  identical, differences on the order of 1e-15) and unconditionally
  O(n).

Combined effect:
[`detect_qrs()`](https://mrpheus.circadia-lab.uk/reference/detect_qrs.md)
on the 468,700-sample recording above now completes in 0.97 seconds
(1855 R-peaks detected, mean HR 118 bpm).

### ✨ New function: `compute_hrv_freq()`

- [`compute_hrv_freq()`](https://mrpheus.circadia-lab.uk/reference/compute_hrv_freq.md)
  — frequency-domain HRV (HF power) locked to each infant’s own
  respiratory rate, rather than the fixed adult-derived HF band
  (0.15–0.4 Hz) that doesn’t apply to neonatal breathing rates (~40–70
  breaths/min at rest). Detects the dominant respiratory frequency
  directly from the physlog’s `resp` channel, resamples the RR-interval
  series via cubic spline, and integrates Welch PSD power in a band
  centered on that detected frequency. Verified against synthetic data:
  exact frequency recovery and correct amplitude-scaling of integrated
  power.
- Fills the gap left by
  [`compute_hrv_sleep()`](https://mrpheus.circadia-lab.uk/reference/compute_hrv_sleep.md),
  which remains a stub expecting staged PSG input (wrong shape for
  physlog data) with its RR-to-HRV step never wired in —
  [`compute_hrv_freq()`](https://mrpheus.circadia-lab.uk/reference/compute_hrv_freq.md)
  is a separate, physlog-shaped path rather than a completion of that
  stub.

### 🧪 Tests

- `test-find-peaks-min-dist.R`: `.find_peaks_min_dist()` against a
  brute-force reference on random and periodic synthetic signals, the
  minimum-distance constraint, tallest-peak-wins tie-breaking (including
  the all-tied flat-signal edge case), and a timing regression guard on
  a large, realistically-smoothed candidate set.
- `test-boxcar-integration.R`: the cumulative-sum boxcar formula against
  `stats::convolve(type = "open")` across several window sizes and
  signal lengths, plus an end-to-end
  [`detect_qrs()`](https://mrpheus.circadia-lab.uk/reference/detect_qrs.md)
  timing smoke test on a synthetic two-minute ECG-like signal.
- `test-compute_hrv_freq.R`: synthetic QRS/respiration generators
  injecting a known respiratory frequency via RSA-modulated RR intervals
  — frequency recovery at two different frequencies, HF power scaling
  with RSA amplitude, band-bounds correctness, the too-few-beats `NA`
  path, and implausible-RR-interval filtering.

All fixtures are fully synthetic — no participant-derived data.

## mrpheus 0.1.4 (2026-07)

### 🚀 Rcpp hot-path optimisation (staging pipeline + event detection)

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

### 🔧 Other changes

- Fixed a latent bug where the pure-R `vapply` definition of
  `.roll_triang_mean` in the normalisation section was shadowing
  `roll_triang_mean_cpp` on every `load_all()`. Now a single Rcpp-backed
  definition exists.
- Removed `zoo` from Imports — no remaining `zoo::` calls anywhere in
  the package.

### 🧪 Tests

All 371 tests pass. `0 errors | 0 warnings | 0 notes`.

## mrpheus 0.1.3 (2026-07)

### 🩺 PSG preprocessing pipeline

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
