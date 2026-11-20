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
  a 3/8-rule variant; scipy (1.13.1 to 1.18.1 checked, see Caveats) uses
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

*Caveats.* YASA 0.8.1 requires scipy >= 1.15 (its PyPI metadata). Measured with
YASA 0.8.1 on a synthetic 2 h recording (240 epochs) under scipy 1.13.1, 1.14.1,
1.15.0, 1.16.3, 1.17.1 and 1.18.1 (numpy 2.5, pandas 2.3.3; 1.13.1 and 1.14.1
are below YASA's own requirement and were run anyway): features and class
probabilities are identical on 1.13.1 to 1.17.1; 1.18.1 differs in one zero
crossing count in the first epoch (EMG channel: 577 against 578; the EOG
channel also changes), which moves one probability by up to 1.1e-3 and no stage; band power
agrees to 1e-13 and the events of `sw_detect` and `spindles_detect` to 6e-14.
Versions before 1.13 were not tested (YASA 0.8.1 needs numpy >= 2.2.2, which
scipy 1.12 and older do not support), so nothing is claimed for them.
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
* `relative` now defaults to `TRUE`, as YASA's does (it was `FALSE`); use
  `relative = FALSE` for absolute power.
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
and `data-raw/generate_bandpower_reference.py`, which used the old settings, and
their reference and report files (`yasa_reference_bandpower.csv`,
`parity_report_bandpower.csv`, `parity_bandpower.png`) are removed.

### 🌊 `compute_slow_oscillations()`: now reproduces `yasa.sw_detect()` (results change)

The detector is a port of `yasa.sw_detect()` (YASA 0.8.1). On the Sleep-EDF
night `SC4001E0` (Fpz-Cz; 8.8 h at 100 Hz, and a 6.5 h version at 256 Hz), with a
stage mask that cuts through events, it finds the same events as YASA: the same
numbers (2,236 and 1,000), event times identical to the sample, and amplitudes
within about 1e-12 uV. The previous version, which paired zero crossings, found
about 20% fewer events than YASA on one night (4.7 against 5.9 per NREM minute,
from a Python mirror of it).

* **Filter.** An FIR band-pass equivalent to `mne.filter.filter_data()` (Hamming
  window, zero phase, 0.3-1.5 Hz with 0.2 Hz transition bands) is applied to the
  whole continuous signal, not epoch by epoch, instead of a Butterworth cascade.
  The design equals `mne.filter.create_filter` to about 1e-16 (checked at 100,
  125, 128, 200, 250, 256, 500 and 512 Hz) and the filtered signal equals
  `filter_data` to about
  5e-13 uV. This also removes the old failure at 256 Hz ("system is
  computationally singular") and the sample-rate limit of the cascade.
* **Detection.** Peaks of the filtered signal are found as `scipy.signal.find_peaks`
  does (negative peaks with a trough of 40-200 uV, positive peaks of 10-150 uV)
  and kept if they lie in the accepted epochs. Each negative peak is paired with
  the closest following positive peak, the peak-to-peak amplitude must be 75-350
  uV, and the durations of the two phases, taken from the zero crossings around
  the peaks, must be 0.3-1.5 s and 0.1-1.0 s. Events that share their start, and
  then their end, with another event are dropped, as in YASA (so two negative
  peaks paired with the same positive peak are both lost). The zero-crossing
  pairing of the old version, and its polarity bug (an epoch whose filtered
  signal started negative returned no events), are gone;
  `detect_so_candidates_cpp()` and its tests are removed.
* **`stages`** are the epochs in which the peaks must lie. As in YASA, an event
  can start or end outside them.
* **Output.** `epoch` is the epoch of the negative peak, and `start_s` and
  `end_s`, relative to it, can lie outside 0-30 s for an event that crosses an
  epoch boundary. New columns: `neg_peak_s`, `mid_crossing_s`, `pos_peak_s`
  (seconds from the start of the first epoch) and `slope_uv_s`.
* **Arguments.** `transition_hz` is now the width of both transition bands of the
  FIR filter (default 0.2, as YASA); zero is rejected, so the old behaviour can no
  longer be requested through it. The other defaults are YASA's, as before.
* A channel whose trimmed standard deviation (central 90% of samples) is outside
  0.1-1000 uV is skipped with a warning, as in YASA.
* Not implemented: YASA's slow oscillation - spindle coupling and
  `remove_outliers`.
* **Tests.** `tests/testthat/test-yasa-slow-oscillation-parity.R` compares the
  result with stored `yasa.sw_detect()` output on a synthetic EEG at 100 and
  256 Hz, the filter taps with `mne.filter.create_filter`, and the detection with
  YASA's own code on a crafted signal that reaches its edge cases (plateaus,
  double-humped waves, limits exactly on the boundary, exact zeros)
  (`data-raw/generate_yasa_slow_oscillation_parity_fixture.py`). Putting each
  deviation back one at a time, 11 of 13 make the tests fail; the other two
  change nothing (removing only the shared-end check, and the wrap-around of the
  previous zero crossing before the first one).

### 💤 `compute_spindles()`: now reproduces `yasa.spindles_detect()` (results change)

The detector is a port of `yasa.spindles_detect()` (YASA 0.8.1). On the Sleep-EDF
night `SC4001E0` (Fpz-Cz; 8.8 h at 100 Hz, and a 6.5 h version at 256 Hz), with a
stage mask that cuts through events, it finds the same spindles as YASA: the same
numbers (73 and 78), start and end times identical to the sample, and every
property within about 1e-11. The previous version (a 11-16 Hz Butterworth
band-pass and an RMS threshold of 1.5 standard deviations, both per epoch, and a
Welch peak frequency) was a different method and had not been compared with YASA.

* **Algorithm.** Two MNE-equivalent FIR filters are applied to the whole
  continuous signal: a broadband (1-30 Hz) and a sigma (12-15 Hz, 1.5 Hz
  transition bands). Three criteria are evaluated at every sample: the relative
  sigma power of a short-time Fourier transform (2 s windows every 0.2 s), the
  moving correlation between the sigma and the broadband signal (0.3 s windows)
  and the moving RMS of the sigma signal (threshold: mean plus 1.5 trimmed
  standard deviations within the accepted epochs, at most 10 uV). The number of
  criteria met is smoothed over 0.1 s, events where all are met are kept, events
  closer than 500 ms are merged, and durations of 0.5-2 s are accepted.
* **Arguments.** `freq_spindle` now defaults to `c(12, 15)` (was `c(11, 16)`).
  New: `freq_broad`, `duration_s` (replaces `min_duration_s` and `max_duration_s`;
  the maximum is 2 s, was 3 s), `min_distance_ms`, `thresh_rel_pow`,
  `thresh_corr` and `thresh_rms` (replaces `threshold_sd`; `NULL` switches a
  criterion off). `rms_window_s` is removed: the window is fixed at 0.3 s, as in
  YASA.
* **Output.** `epoch` is the epoch of the first sample and `start_s` / `end_s`
  are relative to it (`end_s` can exceed 30 s). `peak_freq_hz` is now the median
  instantaneous frequency (was the Welch peak) and `rms_uv` the RMS of the
  detrended broadband signal (was the mean RMS envelope). New columns: `peak_s`,
  `amplitude_uv`, `amp_filtered_uv`, `abs_power`, `rel_power`, `oscillations`,
  `symmetry`.
* A channel whose trimmed standard deviation (central 90% of samples) is outside
  0.1-1000 uV is skipped with a warning, as in YASA.
* Not implemented: YASA's `multi_only`, `remove_outliers` and spindle - slow
  oscillation coupling. `roll_rms_cpp()` and its tests are removed.
* **FIR filter length.** `.mne_fir_bandpass_design()` rounds the filter length up
  to whole samples, as `mne.filter` does; it had rounded to nearest, which gives
  the same length except at rates where the length is not a whole number (for
  example 125 Hz, 2063 against 2062 taps for the slow oscillation filter).
  Checked against `mne.filter.create_filter` for the slow oscillation, spindle
  broadband and sigma filters at 100, 125, 128, 200, 250, 256, 500 and 512 Hz.
* **Tests.** `tests/testthat/test-yasa-spindle-parity.R` compares the result with
  stored `yasa.spindles_detect()` output on a synthetic EEG with spindles at 100
  and 256 Hz, and at 100 Hz on two variants (strong sigma activity in an excluded
  epoch; the signal times 4), plus YASA's intermediate features (STFT relative
  power, moving correlation and RMS)
  (`data-raw/generate_yasa_spindle_parity_fixture.py`);
  `test-compute_spindles.R` tests the helpers against scipy and YASA values.
  Putting each of 19 deviations back one at a time makes the tests fail in every
  case; the cap of 10 uV on the RMS threshold is caught only by its direct unit
  test.

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
