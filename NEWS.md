## mrpheus 0.1.4.9000  (development)

### 🛏️ `stage_epochs()`: now reproduces YASA (results change)

`stage_epochs()` is checked against `yasa.SleepStaging` (YASA 0.8.1, scipy
1.18.1) feature by feature and probability by probability. Three errors in the
feature code, a missing model selection and a missing precision step are fixed;
on the public night used to check this, the stage is now the same as YASA's in
every epoch.

* **Channel-specific classifiers.** YASA picks one of eight classifiers from the
  channels given (EEG, EEG + EOG, EEG + EMG, EEG + EOG + EMG, each with or
  without demographics). mrpheus shipped only the EEG + EOG + EMG model and
  filled missing channels with NA. With EOG absent that gave 50% accuracy
  against manual scoring on the test night (YASA: 83%); with EMG absent it
  differed from YASA's EEG + EOG model in 91 of 780 epochs. All eight are now
  bundled as `inst/models/yasa_<channels>[_demo].txt.gz` (about 3.3 MB each) and
  chosen automatically. New arguments `age` and `male` select the `+demo`
  classifiers; `model_path` now defaults to `NULL` (automatic) and accepts
  `.gz` files. The result carries the classifier as `attr(, "classifier")`.
  Run `data-raw/fetch_yasa_models.py` to create the model files; the old
  `inst/models/yasa_staging.txt` is byte-identical to `yasa_eeg_eog_emg` and
  can be deleted.
* **Single precision.** YASA rounds all features to float32 before predicting;
  mrpheus now does too.
* **Simpson integration for an even number of points.** `.scipy_simpson()` used
  a 3/8-rule variant; scipy (>= 1.11, checked on 1.14.1 and 1.18.1) uses
  Cartwright's correction on the last interval. This changed the 0.4-1 Hz and
  1-4 Hz band powers (even point counts at 0.2 Hz resolution) and the
  delta ratios built from them, by up to half a standard deviation, and also
  affects `compute_band_power()` and `compute_temporal_bandpower()` wherever a
  band holds an even number of frequency bins. The old tests used linear data,
  which every Simpson variant integrates exactly, so they could not see it.
* **Skewness and kurtosis** standardised with the sample (n - 1) standard
  deviation; YASA uses population central moments. The normalised variants
  cancel the difference, the raw ones (which the classifiers also use) did not.
  A flat epoch now gives `NaN` skew and kurtosis, as in YASA.
* **Zero crossings** now follow `antropy.num_zerocross` (`np.signbit`): an exact
  zero is no longer counted as two crossings.

*Measured agreement.* On Sleep-EDF `SC4001E0` (6.5 h, 780 epochs), 100 Hz and
256 Hz inputs, all 16 classifier configurations give the same stage in every
epoch as YASA, with probabilities equal to about 1e-16 in 99.9% of epochs. The
exception is the recording's first and last sample: the zero-phase FIR filter's
output there is zero by construction and its sign is floating-point noise
(+6e-15 in MNE, -4e-14 here), so one zero-crossing count can differ by one.
Because `robust_scale` uses whole-column quantiles this shifts every epoch's
normalised value slightly; the largest probability difference seen was 1e-3,
in one epoch of 12,480. Taking only the zero-crossing features from YASA makes
the probabilities exactly identical (difference 0.0). `tests/testthat/test-yasa-parity.R`
checks all eight classifiers against a stored YASA reference.
The earlier "validated" parity script compared at MAE 0.01 and skipped 15 edge
epochs; that was not parity to float precision.

*Caveats.* Agreement is with YASA 0.8.1 on scipy >= 1.15; older scipy versions
integrate even-length bands differently and give different YASA features.
Preprocessing with `preprocess_psg()` before staging is not what YASA does and
changes the result (84 of 780 epochs on the test night, mostly through the EOG
band-pass); stage from unpreprocessed signals to reproduce YASA.

### 📊 `compute_band_power()`: now follows `yasa.bandpower()` (results change)

The defaults did not match YASA's, although the help page said they did. On
Sleep-EDF `SC4001E0` (780 epochs, two channels, at 100 and 256 Hz; YASA 0.8.1)
the old defaults differed from `yasa.bandpower()` by a median of 6 to 12% per
band and epoch in absolute power (and in relative power), with a few epochs far
larger. The defaults are now YASA's, and the result agrees with
`yasa.bandpower()` to about 1e-14 relative in absolute power and about 1e-12 in
relative power (maximum over all epochs, channels and both sampling rates; also
checked against `yasa.bandpower()` itself on ten epochs).

* **Segment averaging:** the bias-corrected median of scipy (`average = "median"`),
  not the mean. This is the largest single difference, about 10% of band power
  at the median and more than 200% at the 99th percentile of epochs.
* **Window:** Hamming, not Hann (about 2%), and the periodic form that
  `scipy.signal.welch` uses. mrpheus built the symmetric form, which alone gave up
  to 0.8% in absolute band power with every other setting matched. New argument
  `window` (`"hamming"` or `"hann"`).
* **`nfft` and `noverlap`** default to `NULL`: the window length (YASA's, scipy's)
  and half a window. The old defaults of 1024 and 200 samples differed from that
  by about 3% (nfft, at 100 Hz) and 7% (overlap, at 256 Hz); 200 samples is half a
  window only at 100 Hz.
* **Gamma** is 30-40 Hz (YASA's default bands), not 30-45 Hz.
* **Relative power** is a band's Simpson integral divided by the PSD integral
  between the lowest and highest band edge, as in YASA, and no longer divided by
  the sum of the band powers. The two are equal when the bands tile that range
  and their edges fall on FFT bins (the new defaults), but differed by a median
  of 0.7% and up to 5.7% at 100 Hz with `nfft = 1024`. `total_power` is now that
  integral (YASA's `TotalAbsPow`), not the sum of the band powers.
* `relative` still defaults to `FALSE`; YASA's default is `TRUE`.
* The old estimate can be requested with `window = "hann", average = "mean",
  noverlap = 200L, nfft = 1024L` and `bands` with gamma `c(30, 45)`, except that
  the window is now the periodic form.

`compute_temporal_bandpower()` is unchanged (its output is identical to before):
it has no YASA counterpart, it keeps its own settings (symmetric Hann window,
mean, gamma 30-45 Hz) and agrees with an independent scipy implementation of
the same method to about 1e-13. Its help page no longer says that its default
bands match `compute_band_power()`.

`tests/testthat/test-yasa-bandpower-parity.R` compares `compute_band_power()`
with a stored `yasa.bandpower()` reference (`data-raw/generate_yasa_bandpower_parity_fixture.py`);
each of the old settings above makes it fail, apart from the relative-power
definition, which `test-compute_band_power.R` checks with bands that leave a gap.

The "Parity with YASA" section of the spectral-analysis vignette and the README
wording are rewritten. They described the old settings as validated, and put the
~0.75% relative-power gap down to floating-point accumulation; it was the
relative-power definition. The vignette now reports the measured agreement on
SC4001E0 (2,650 epochs): at most 8e-15 relative in absolute power and 1e-12 in
relative power. The old validation scripts `data-raw/validate_band_power_parity.R`
and `data-raw/generate_bandpower_reference.py`, which used the old settings, are
removed.

### ⚠️ `compute_slow_oscillations()`: fixes and new defaults (results change)

* **Polarity bug fixed.** `detect_so_candidates_cpp()` paired zero crossings
  in steps of two starting from the first crossing in the epoch, whatever its
  direction. In an epoch whose filtered signal started negative every pair was
  half a cycle out of phase and the epoch returned no slow oscillations, so
  roughly half the epochs could be lost. The walk now starts at a crossing
  after which the signal is negative. Existing tests passed because their
  signals started at exactly 0.
* **Filter fixed at higher sample rates.** The single 8th-order band-pass was
  numerically singular at 256 Hz (`gsignal::filter_zi()` failed with
  "system is computationally singular"). It is now a zero-phase high-pass /
  low-pass cascade, designed once rather than per epoch.
* **New defaults, taken from `yasa.sw_detect()`:** `freq_so = c(0.3, 1.5)`
  (was `c(0.5, 2)`), `duration_neg_s = c(0.3, 1.5)` (was `c(0.1, 1.5)`),
  `amp_ptp_threshold_uv = c(75, 350)` (was `c(75, 500)`). New arguments
  `amp_neg_uv = c(40, 200)` and `amp_pos_uv = c(10, 150)` (trough and peak
  limits; `NULL` switches them off) and `transition_hz = 0.2`, which places
  the filter's -6 dB points at 0.2 and 1.6 Hz as YASA's FIR does. The old
  behaviour can be requested with `freq_so = c(0.5, 2)`,
  `duration_neg_s = c(0.1, 1.5)`, `amp_ptp_threshold_uv = c(75, 500)`,
  `amp_neg_uv = NULL`, `amp_pos_uv = NULL`, `transition_hz = 0`.
* **Not identical to YASA.** The pairing is still zero-crossing based, and
  epochs are filtered one at a time. In a Python comparison on one Sleep-EDF
  night (SC4001E0, Fpz-Cz, N2 + N3), a mirror of the old defaults gave 8.7
  events per NREM minute against 5.9 for YASA, with 53% of its events
  matching a YASA event (negative peaks within 0.5 s). A mirror of the new
  defaults gave 4.7 per minute, with 90% matching a YASA event and 72% of
  YASA's events found. This is one night, one channel and a Python mirror,
  not a validation.

### 🔧 Resampler and spectrogram fixes

* **`.mne_fft_resample()` now downsamples.** Staging resamples every channel to
  100 Hz, so 256 Hz input (a common PSG rate) stopped with `only 0's may be
  mixed with negative subscripts`: the function only implemented the upsampling
  path, and the Sleep-EDF data used so far is already at 100 Hz. For a shorter
  output it now truncates the spectrum as MNE's `_fft_resample` does (Nyquist
  bin doubled for an even new length, bins scaled by `new_len / orig_len`).
  A Python mirror of the new logic differs from `mne.filter.resample` (MNE
  1.13.2) by about 5e-14 on four lengths (7,680 to 460,800 samples). The
  upsampling code is unchanged. Tests: `tests/testthat/test-mne-fft-resample.R`.
* **`compute_spectrogram()` overlap.** The function passed `noverlap` to
  `gsignal::pwelch()`, which has no such argument (it takes `overlap`, a
  fraction), so the call stopped with "unused argument". It now passes its own
  `overlap` argument (default `0.5`) straight through. The function had no
  tests; `tests/testthat/test-compute_spectrogram.R` adds a synthetic 10 Hz
  check (it runs, returns one row per epoch, and peaks within 0.5 Hz of 10 Hz).

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
